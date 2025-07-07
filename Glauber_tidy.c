#include <assert.h>
#include <stdlib.h>
//#include <omp.h>
#include <time.h>
#include <string.h>   // for memset
#include<math.h>
#include<gsl/gsl_rng.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <ctype.h>


#define MAX_LINE 256


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////// This code simulates the continuous time Glauber Dynamics //////////////////////////////// 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The update is sinchronous. This leads to a problem
// At an interface, if both the right and the left spin flip, then two interfaces are created. This seems
// not to lead to any problem if the interfaces are anihilated at a rate of order 1 at each microscopic step 
// If the annihilation rate is 1/N^2, then this creation rates takes the lead and the invariant meausre  is not the 
// Ising one. To solve this problem, here the Glauber update is performed in only at odd or evens sites.

// Please Observe that for the discrete time Glauber dynamics N^2*prob_creation = prob_annihilation 
// while here prob_annihilation = 1 - exp(annihilation_rate/n_steps) prob_creation = 1- exp(creation_rate/n_steps)
// with annihilation_rate = creation_rate *N^2, resulting in a slightly smaller annihilation probability with respect to the 
//discrete Glauber Dynamics

// STRUCTURES 

typedef void (*Initialization_function)(int *spin, int l, void *variables);  // This pointer to a function will point to the function that initializes the spins at time 0. 
// The void *variable pointer will point to the variable that the initialization function needs, and it is void since the set of variables is different for each function.

/// CONFIG CONTAINS ALL THE PARAMETERS OF THE SYSTEM, INCLUDING THE INITIALIZATION FUNCTION 
typedef struct { 
	// The parameters needed to define
	int N; // scaling parameter
	double T;
	int L; 
	double annihilation;  // annihilation rates of the microscopic model
	double creation;  // creation rates of the microscopical model
	// The flipping rates for an interface spin are assumed to be 1; 
	//  The parameters needed to discretize the dynamics

	int N_simulations; // number of independent simulations
	int resolution; // Number of subntervals in which to divide the unit of macroscopic time. Since we need N^2 time steps of the microscopic dynamics, we need to plot N*N/resolution micro evolutions
	int Micro_n_steps; // Number of time steps in which to divide the microscopic dynamics. It should be made in such a way that annihilation/Micro_Time_Step << 1, and we must pay attention to the resolution for the Mersenne Twister, that is creation/Micro_Time_Step must be detectable by a uniformly generated random variable.  
	
	Initialization_function initialize; // Pointer to the function that initializes the spins at time 0.
	void *initialization_param; // The void pointer will hold the third argument of the above function. It will be a double pointer if the initial distribution needs a parameter, a NULL pointer otherwise
	// Parameters needed to compute the statistics 
	
	//Average Magnetization in an interval 
	float A1_magnetization; 
	float B1_magnetization; 
	// Correlation Estimation 
	float A1_corr; 
	float A2_corr; 
	float B1_corr; 
	float B2_corr; 
	// Time Delayed Correlation Estimation
	float A1_time_delayed_corr; 
	float A2_time_delayed_corr; 
	float B1_time_delayed_corr; 
	float B2_time_delayed_corr; 
//	float Lambda_T; 
	} Configuration; 

/// THE FOLLOWING STRUCTURES CONTAIN THE PARAMETERS TO BE PASSED TO DIFFERENT FUNCTIONS 

// To be passed to the Glauber update function

typedef struct {
	int l; 
	double prob_annihilation;  // 1- exp( - rate * dt)
	double prob_creation; 
	double prob_movement; 
	int n_micro_evolution; // Number of times we evolve the dynamics before reaching the required resolution (N*N*Micro_n_steps)
	gsl_rng *rng;
	}  Parameters_update ;

// Containing the intervals where to compute the correlations and the magnetizations 

typedef struct {
	int a1; //left endpoint of the first microscopic interval
	int b1; //left endpoint of the second microscopic interval 
//	int l; // Theoretically, this parameters should go in the definition of spins, which contains also the neighbour topology.
	}  Interval;  

typedef struct {
	int a1; //left endpoint of the first microscopic interval
	int b1; //left endpoint of the second microscopic interval 
	int l; // Theoretically, this parameters should go in the definition of spins, which contains also the neighbour topology.
	}  Interval_with_l;  

typedef struct {
	int a1; //left endpoint of the first microscopic interval
	int a2; //left endpoint of the second microscopic interval 
	int b1; // right endpoint of the first microscopic interval 
	int b2; // right endpoint of the second microscopic interval 
//	int l; // Theoretically, this parameters should go in the definition of spins, which contains also the neighbour topology.
	}  Two_Intervals; 


// FUNNCTIONS BETWEEN STRUCTURES 

// Function that set the intervals whete to compute the magnetization, n_interfaces, time_delayed correlations and correlations

void set_parameters( Configuration  *conf, Parameters_update *parameters );

void set_interval_magnetization(Configuration *conf, Interval *interval_extrema){
	interval_extrema -> a1 = (int) ((conf -> N)*(conf -> A1_magnetization)); 
	interval_extrema -> b1 = (int)( (conf -> N)*(conf ->B1_magnetization )); 
//	interval_extrema -> l = (int)( conf -> N * conf -> L ); 
}

//void set_interval_n_interfaces(Configuration *conf, Interval_with_l *interval_extrema){
//	interval_extrema -> a1 = (int) (conf -> N)*(conf -> A1_n_interfaces); 
//	interval_extrema -> b1 = (int) (conf -> N)*(conf ->B1_n_interfaces); 
//	interval_extrema -> l = (int )( conf -> N * conf -> L );
//}

void set_time_delayed_params(Configuration *conf, Two_Intervals *Interval_extrema){
	Interval_extrema -> a1  = (int) ( (conf -> N)*(conf -> A1_time_delayed_corr));  
	Interval_extrema -> a2  = (int) ((conf -> N)*(conf -> A2_time_delayed_corr));  
	Interval_extrema -> b1  = (int) ((conf -> N)*(conf -> B1_time_delayed_corr));  
	Interval_extrema -> b2  = (int) ((conf -> N)*(conf -> B2_time_delayed_corr)); 
//	Interval_extrema -> l = (int)( conf -> N * conf -> L );  	
}
void set_correlation_params(Configuration *conf, Two_Intervals *Interval_extrema ) { // Some checks when reading the conf file REMEMBER 
	Interval_extrema  -> a1  =(int) ((conf -> N)*(conf -> A1_corr));
	Interval_extrema  -> a2  =(int) ((conf -> N)*(conf -> A2_corr));  
	Interval_extrema  -> b1  =(int) ((conf -> N)*(conf -> B1_corr));  
	Interval_extrema  -> b2  =(int) ((conf -> N)*(conf -> B2_corr)); 
//	Interval_extrema -> l = (int ) (conf -> N )* (conf -> L);  
}


void write_stats_header(FILE *f, int N_sim, const char **stat_names, int N_stats) {
    fprintf(f, "time");
    for (int i = 0; i < N_sim; i++) {
        for (int j = 0; j < N_stats; j++) {
            fprintf(f, "\tsim%d_%s", i, stat_names[j]);
        }
    }
    fprintf(f, "\n");
}
void write_stats_row(FILE *f, double *corr, double *magnetization, double *time_delayed_correlations, int *interfaces, int N_sim) {
    for (int i = 0; i < N_sim; i++) {
            fprintf(f, "%.6f\t", corr[i]);
            fprintf(f, "%.6f\t", magnetization[i]);
            fprintf(f, "%.6f\t", time_delayed_correlations[i]);
	    fprintf(f, "%d\t", interfaces[i]);	
	 }
    fprintf(f, "\n");
}
// FUNCTIONS AND STRUCTURES FOR THE STRUCTURE CONFIGURATION 
// A MESS TO READ THE CONFIGURATION FROM A FILE

typedef int (*ConfigSetter)(Configuration *conf, const char *value); // The function that associate the character "5"  conf -> some_field_of_conf = 5 will be of this type. They will return the successfullness of the operation 



int set_N(Configuration  *conf, const char *value);
int set_T(Configuration  *conf, const char *value);
int set_L(Configuration  *conf, const char *value);
int set_N_simulations(Configuration *conf, const char *value); 
int set_resolution(Configuration *conf, const char *value);
int set_Micro_n_steps(Configuration *conf, const char *value);
int set_A1_magnetization(Configuration *conf, const char *value);
int set_B1_magnetization(Configuration *conf, const char *value);
int set_A1_corr(Configuration  *conf, const char *value); 
int set_A2_corr(Configuration  *conf, const char *value); 
int set_B1_corr(Configuration  *conf, const char *value); 
int set_B2_corr(Configuration  *conf, const char *value);
int set_A1_time_delayed_corr(Configuration *conf, const char *value); 
int set_A2_time_delayed_corr(Configuration *conf, const char *value);
int set_B1_time_delayed_corr(Configuration *conf, const char *value); 
int set_B2_time_delayed_corr(Configuration *conf, const char *value);
//int set_Lambda_T(Configuration *conf, const char *value);


typedef int (*ConfigSetter)(Configuration *conf, const char *value); // The function that associate the character "5"  conf -> some_field_of_conf = 5 will be of this type. They will return the successfullness of the operation 
typedef struct {
    const char *key;
    ConfigSetter setter;
} ConfigEntry; // This structure will be used to associate to each character which field of conf to fill, that is, which  Config_setter to use 

ConfigEntry config_table[] = {
    {"N", set_N},
    {"T", set_T},
    {"L", set_L},
    {"N_simulations", set_N_simulations}, 
    {"resolution", set_resolution},
    {"Micro_n_steps", set_Micro_n_steps},
    {"A1_magnetization", set_A1_magnetization},
    {"B1_magnetization", set_B1_magnetization},
    {"A1_corr", set_A1_corr},
    {"A2_corr", set_A2_corr},
    {"B1_corr", set_B1_corr},
    {"B2_corr", set_B2_corr},
    {"A1_time_delayed_corr", set_A1_time_delayed_corr},
    {"A2_time_delayed_corr", set_A2_time_delayed_corr},
    {"B1_time_delayed_corr", set_B1_time_delayed_corr},
    {"B2_time_delayed_corr", set_B2_time_delayed_corr},
//  {"Lambda_T", set_Lambda_T},
    {NULL, NULL} //A sentinel indicating the end of the structure. 
}; // I need to validate the variables that I read!!!REMEMBER Ok, let's use the simple suggestion  and write a validating function 

int is_scaled_integer(float val, int N) {
    return fabsf(val * N - roundf(val * N)) < 1e-5;
}

void validate_config(const Configuration *conf);

void read_config_file(Configuration *conf, const char *filename);

void read_config(Configuration *conf ); //Reads the configuration parameters that define the dynamics. It calls the following two functions  
	
void choose_rates(double *annihilation, double *creation, int N );// Function to read the rates used in the simulations  

void choose_initialization_function(Configuration *conf); // Function that reads the function that will be used to initialize the spins, between the initialize spin functions defined later 

void print_config_to_screen(Configuration *configuration); // Prints the chosen parameters to the keyboard, with the exception of the initial distribution 

void print_conf_json(const Configuration *conf, const char *dirname);


// FUNCTIONS FOR THE TYPE PARAMETERS 


void free_parameters_update(Parameters_update *par){
	free( par -> rng); 
	free( par); 
} //frees the memory


// FUNCTIONS FOR THE OUTPUT SYSTEM CALL  
// A function needed since fgets in the function choose rates gave me odd results 
void clear_stdin() {
    int c;
    while ((c = getchar()) != '\n' && c != EOF);
}


void create_output_directory(const Configuration  *conf, char *dirname_out, size_t len);


int is_comment_or_blank(const char *line);



// We initialize the spins at time zero using one of the following functions 
void initialize_spins_1(int *spins, int l, void *params); // Geometric jumps between spin +1 and -1  of parameter lambda/N, lambda to be chosen. If Lambda = 1, it is the invariant measure for the Glauber rates 
void initialize_spins_2(int *spins, int l, void *params); // all -1 
void initialize_spins_3(int *spins, int l, void *params); // Two interfaces
void initialize_spins_4(int *spins, int l, void *params); // Each spin drawn indepednently from the other ones, +1 with probability p to be chosen
void initialize_spins_5(int *spins, int l,  void *params );
//void initialize_spins_6(int *spins, int l, void *params); 
//void initialize_spins_7(int *spins, int l, int k, double p, int tot_sites); // Invariant measure, as in 1, coditioned on having spins +1 on [L/2, L/2+k]

//void Glauber_update(int *spins_aux , int *new_spins, Parameters_update *params);
 
void Glauber_update_parity(int *spins_aux , int *new_spins, Parameters_update *params, int a );

double compute_magnetization(int *spin, Interval *interval);
void compute_and_print_n_interfaces (FILE *f, int data, int tot_sites, int l   ); 
int n_interfaces(int *spins, int l ); 
//void print_spins_to_file_binary(int *spins,  FILE *f, Parameters_update *param); 

void print_spins_to_file_binary2(int *spins,  FILE *f, int length); 

double F(double x, int N); 

double GF( int l, int N, int *spins );

int creation_of_interfaces(int a1, int b1, int l,  int *spin1, int *spin2);

double correlation( int *spin1, int *spin2, Two_Intervals *intervals); 




//////////////////////////////////////////////////////////////////////////////////////
// MAIN 
//////////////////////////////////////////////////////////////////////////////////////


int main (){
// Allocate the variables defining the Dynamics and the Statistics to be computed

	
	Configuration *conf; 	
	conf = malloc(sizeof(Configuration)); 
        	if (conf == NULL) {
      			 perror("conf malloc failed");
   	       		  exit(EXIT_FAILURE);
   		 }
		
	// Read from the keyboard the variables  of the program	
	read_config_file(conf, "config.txt"); 	
	print_config_to_screen(conf);	

	///////// CREATE THE OUTPUT DIRECTORIES AND FILES 
	// Variable that saves the output directory 
	char dirname[256];  // You can adjust size as needed
	// Create the output directory 
	create_output_directory(conf, dirname, sizeof(dirname)); 	
	// Print the choices in a file json
	print_conf_json(conf, dirname);
	int timestep_marker = -9999; // Marker for the binary output file 

	// Defining the parameters to the variables used to compute the functions 
	// Allocate the parameters to compute statistic 
	Parameters_update *parameters; 	
	Two_Intervals *corr_par; 
	Two_Intervals *time_delayed_corr_par; 	
	Interval *par_av_magnetization; 
	
	 parameters = malloc(sizeof(Parameters_update));  
	 if ( parameters == NULL) {
      		 perror("parameters malloc failed");
   	         exit(EXIT_FAILURE);
   	 }
	corr_par = malloc(sizeof(Two_Intervals));  
        if (corr_par == NULL) {
      		 perror("Correlation Parameters malloc failed");
   	         exit(EXIT_FAILURE);
   	 }
	time_delayed_corr_par = malloc(sizeof(Two_Intervals));  
        if (time_delayed_corr_par == NULL) {
      		 perror("Time Delayed Correlation Malloc  failed");
   	         exit(EXIT_FAILURE);
   	 }
	par_av_magnetization = malloc(sizeof(Interval));  
        if (par_av_magnetization == NULL) {
      		 perror(" Average magnetization malloc failed");
   	         exit(EXIT_FAILURE);
   	 }
	// SETTING THE VALUES

	// Set the parameters to be passed to the update function of Glauber dynamics
	set_parameters(conf, parameters);
   	set_interval_magnetization(conf, par_av_magnetization); 
	set_correlation_params(conf, corr_par);
	set_time_delayed_params(conf, time_delayed_corr_par); 

	// Allocating the variables where to save the statistics
	double *corr = malloc((conf -> N_simulations)*sizeof(double)); // SAFEY CHECK REQUIRED
	double *magnetization = malloc((conf -> N_simulations)*sizeof(double)); // SAFEY CHECK REQUIRED
	int *interfaces = malloc((conf -> N_simulations)*sizeof(double)); // SAFEY CHECK REQUIRED
	double *time_delayed_corr = malloc((conf -> N_simulations)*sizeof(double)); // SAFEY CHECK REQUIRED
//	double *generator = malloc((conf -> N_simulations)*sizeof(double)); // SAFEY CHECK REQUIRED
	
	int length = (conf -> N )*(conf -> L ) * (conf -> N_simulations); 
 
	printf("  At each microscopic time_step, the flipping probabilities are \n" ); 
	printf(" For an interface spin: %lf \n",parameters -> prob_movement);
	printf(" For a parallel spin in the middle of the sea (creation): %lf \n", parameters -> prob_creation);
	printf(" For an antiparallel spin in the middle of the sea (annihilation): %lf \n",  parameters -> prob_annihilation); 	
	printf("DEBUG: N = %d, T = %lf, A1_corr, A2_corr , B1_corr, B2_corr = %d, %d, %d, %d", conf -> N, conf -> T, corr_par -> a1,  corr_par -> a2, corr_par -> b1, corr_par -> b2 );
	
	// ALLOCATE THE SPINS 
	// We allocate the contiguous memory slots for the spins, and we do not use matrices since it is easy to separate it 
	int *data = malloc(length*sizeof(int)); 
	int *new_data = malloc( length *sizeof(int));
	int *past_data = malloc(length*sizeof(int));
	int *data_aux = malloc( length*sizeof(int));
	
	// We now allocate the pointers to the initial spisn 

	int **spins = malloc((conf -> N_simulations)*sizeof(int * )); //SANITYCHECK
	for(int i = 0; i < conf -> N_simulations; i++){
		spins[i] = &data[i *(parameters -> l)]; 
	}

	int **first_spins = malloc((conf -> N_simulations)*sizeof( int * ));//SANITY CHECK
	for(int i = 0; i < conf -> N_simulations; i++){
		first_spins[i] = &past_data[i *(parameters -> l)]; 
	}
	
	int **spins_aux = malloc((conf -> N_simulations)*sizeof(int *));//SANITY CHECK
	for(int i = 0; i < conf -> N_simulations; i++){
		spins_aux[i] = &data_aux[i *(parameters -> l)]; 
	}
	
	int **new_spins = malloc((conf -> N_simulations)*sizeof(int *));// SANITY CHECK REQUIRED 
	for(int i = 0; i < conf -> N_simulations; i++){
		new_spins[i] = &new_data[i *(parameters -> l)]; 
	}

/// INITIALIZE THE SPINS 
	int k = conf -> N;
//	initialize_spins_6(data,  parameters -> l,k , 1.0, length); 
	for(int i = 0; i < (conf -> N_simulations); i++ ){	
	conf -> initialize(spins[i], parameters -> l, conf -> initialization_param); 	
	}	
	memcpy(past_data, data, length*sizeof(int)); // We initialize the past data  
	memcpy(new_data, data, length*sizeof(int)); 
      // Here we decompose the for loops in the following way: The total number of steps is N*N*Time_Steps*T. We decompose them in such a way that we hide n_micro_steps in
      // the update function. When we do an update we plot the microscopic configuration. When we do N*N/(n_micro_steps) of such increments we will have made a macroscopic increment of 1/T. 
      // and we now can compute the one time statistics and the statistics that depend on the infinitesimal increments. 
      // When we do Lambda_T/ Time_step of such increments, we can compute the two time statistic, when the time is at distantce Lambda_T
      //  we can perform this T/Lambda_T times. 

      // What are the relations to be satisfied if I want to write precise for loops? well, N*N/(n_micro_time_steps) is integer, I don't think is necessary that T/Lamda_T is an integer, but maybe it would be nice that Lambda_T/ dT = Time_Steps * Lambda_T is an integer. So this are the conditions. 

	// Open the file in the output folder 
	char filepath[512];
	snprintf(filepath, sizeof(filepath), "%s/spins_output.bin", dirname);
	FILE *output = fopen(filepath, "wb");
	if (!output) {
		perror("File open failed");
		exit(1);
	}
	
	char filepath2[512];
	snprintf(filepath2, sizeof(filepath2), "%s/Statistics.txt", dirname);
	FILE *Statistics = fopen(filepath2, "w");
	if (!Statistics) {
		perror("File open failed");
		exit(1);
	}
	const char *stat_names[] = {"correlations", "time_delayed_corr", "magnetization", "interfaces"};
	int N_stats = 4;
	write_stats_header( Statistics, conf -> N_simulations,  stat_names, N_stats); 


	//HERE WE INITIALIZE THE VARIABLES THAT WILL BE UPDATED 	
	int a = 0; 
	int Macro_Times; 
	Macro_Times =(int) ((conf -> T)*(conf -> resolution));  
	for (int t = 0; t < Macro_Times; t++ ){		
      		fprintf(Statistics, "%lf\t", ((double) t)/((double) conf -> resolution));  
		for(int i = 0;  i<  conf -> N_simulations; i++){
		// GLAUBER UPDATE
			Glauber_update_parity(spins_aux[i], new_spins[i] , parameters, a ); // Here I need the probabilities, the rng, N*L, l, n_evolve = (Micro_n_steps*N*N)/resolution    	
			// Here I compute the statistic 
			corr[i] = correlation(spins[i],spins[i], corr_par); // This is the fixed time correlation. !! With this choice, The last one is not plotted!! 
			magnetization[i] = compute_magnetization(spins[i], par_av_magnetization); // This is the magnetization
			time_delayed_corr[i] = correlation(first_spins[i], new_spins[i], time_delayed_corr_par);
			interfaces[i] = n_interfaces(spins[i], parameters -> l ); 
		}
		print_spins_to_file_binary2(data, output, length); 
		fwrite(&timestep_marker, sizeof(int), 1, output);
		write_stats_row(Statistics, corr, magnetization, time_delayed_corr, interfaces, conf -> N_simulations); 	
		memcpy(data, new_data, length*sizeof(int)); // The data is only needed in order to compute the infinitesimal increments.  
		a = (a + 1) % 2; 
//		fprintf(Generator, " %lf",  ((double) t)/((double) conf -> T)); 
//		for(int i = 0; i < (conf -> N_simulations); i++){
//			fprintf(Generator, " %lf", GF( conf -> parameters -> l, conf -> N, (new_data + (conf -> parameters -> l )*i )) - GF(conf -> parameters -> l, conf -> N,  (data + (conf -> parameters -> l)*i) ) ); 		
//		}
//		fprintf(Generator, "\n");
		

	} 
	fclose(output); 
	fclose(Statistics);
//// CLEANING THE MEMORY 
	
	free(data); 
	free(new_data); 
	free(past_data); 
	free(data_aux); 
	free(spins); 
	free(new_spins); 
	free(spins_aux); 
	free(first_spins); 
	free_parameters_update(parameters);
	free(corr_par);
	free(time_delayed_corr_par);
	free(par_av_magnetization); 
	free(conf); 

	// HERE I CALL THE PYTHON SCRIPT THAT PLOTS THE SPIN HISTORY 
	char command[300];
	snprintf(command, sizeof(command), "python3 Plot_Glauber_tidy_Binary.py \"%s\"", dirname);
	int ret = system(command);
	if (ret != 0) {
		fprintf(stderr, "Failed to run Plot_Glauber_tidy.py on %s\n", dirname);
	}
	// Now I directly open the spin simulations 
	char open_cmd[300];
	snprintf(open_cmd, sizeof(open_cmd), "open %s/si* ", dirname);
	system(open_cmd);
	
	// HERE I CALL THE PYTHON SCRIPT THAT PLOTS ME THE STATISTICS 
	char command2[300];
	snprintf(command2, sizeof(command2), "python3 Plot_stats.py  ./%s/statistics.txt  --outdir ./%s/", dirname, dirname);
	int ret2 = system(command2);
	if (ret != 0) {
		fprintf(stderr, "Failed to run Plot_Stats.py on %s\n",dirname);
	}
	// Now I directly open the spin simulations 
	char open_cmd2[300];
	snprintf(open_cmd2, sizeof(open_cmd2), "open %s/*.png ", dirname);
	system(open_cmd2);
	

	return 0;
}





double correlation( int *spins1, int  *spins2, Two_Intervals *Intervals){
	double corr = 0;
//DEBUG
assert(Intervals != NULL); 
        for(int i = Intervals -> a1; i < Intervals -> b1 ; i++ ){
		for( int j = Intervals -> a2; j < Intervals -> b2; j++)
		corr = corr + spins1[i]*spins2[j]; 
	}
	return corr/((double) ( ((Intervals -> b1) - (Intervals -> a1) )*((Intervals -> b2) - (Intervals ->  a1) )));
}


double compute_magnetization(int *spin, Interval *interval){
	double sum = 0; 	
	for(int i = interval -> a1; i < interval -> b1; i++ ){
		sum = sum + spin[i]; 
	}
	return sum/((double ) ((interval -> b1) - (interval -> a1))); 
}

//int n_interfaces(int *spin, Interval_with_l *interval){
//	int sum = 0; 
//	for(int i = interval -> a1; i < interval -> b1; i++){
//		sum = sum + (1 - spin[i]*spin[i+1])/2;
//	}
//	if( interval -> b1 == interval -> l ){
//		sum = sum +(1 - spin[interval -> l - 1] * spin[0])/2; 
//	} else {
//		sum = sum + (1- spin[interval -> b1]*spin[interval -> b1 + 1])/2; 
//	}
//	return sum; 	
//}

int creation_of_interfaces(int a1, int b1, int l, int *spin1, int *spin2){
	int num_interfaces = 0; 
	for (int i = a1; i < b1; i++){
		if ( spin1[i] != spin1[(i+1 + l) % l ]){
			num_interfaces = num_interfaces + 1;
		}
		if (spin2[i] != spin2[(i+1 +l )% l]) {
		num_interfaces = num_interfaces - 1;
		}
	}
	return (num_interfaces); 
}

int n_interfaces( int *spins, int l){
	int sum = 0; 
	if (spins[0]!=spins[l-1]){
		sum = 1; 	
	}
	for(int i = 0;  i < l-1; i++){
		if( spins[i ] !=spins[i+1]){
			sum = sum +1; 		
		}
	}
	return sum;
}

double F(double x, int N) {
    return sin( 4* M_PI *x / N);  // microscopic coordinate x ∈ [0, L]
} // Check if I need to put some boundary conditions. This function should be an integral of an L^\infty function over the torus 

double GF(int l, int N,  int *spins ){
	double integral = 0;
	for( int i = 0; i < l;  i++){
		integral = integral + spins[i]*F(i/N, N); 
	}
return integral/((double) N ); 	
}

//int count_interfaces(int *spin, int N, double a, double b) {
//    int start = (int)floor(a * N);
//    int end = (int)floor(b * N);
//    int interfaces = 0;
//
//    // Make sure we correctly wrap around the torus
//    for (int i = start; i < end; ++i) {
//        int idx1 = i % N;
//        int idx2 = (i + 1) % N;
//        if (spin[idx1] != spin[idx2]) {
//            interfaces++;
//        }
//    }
//
//    return interfaces;
//}

//
//void compute_and_print_n_interfaces(FILE *f,int  data, int N_simulations, int l, int a1, int b1 ){
//
//	for(int  j = 0; j < N_simulations; j++ ){
//		int sum = 0; 
//		for(int i = a1; i < b1; i++){
//			sum = sum + data[j*l + i ]; 
//		}
//		fprintf(f, "%lf", (double) sum/ ((double) params -> N))
//	}
////%} 


void read_config(Configuration *conf ){
	printf(" Choose the scaling parameter N \n"); 
	scanf("%d",  &conf -> N ); 
	printf(" Choose the final macroscopic time T \n"); 
	scanf("%lf", &conf -> T ); 
	printf(" Choose the macroscopic size of the torus L \n" ); 
	scanf("%d", &conf -> L); 
	// The time step 
	printf(" Choose the macroscopic time resolution (number of plotted times )/(Macroscopic time) \n");
	scanf("%d", &conf -> resolution);
	printf(" Choose the number of intervals n which to divide the microscopic unit of time \n"); // It must be that resolution = N*N *Micro_time_step
	int valid = 0; 
	do {
		scanf("%d", &conf -> Micro_n_steps);
		if( (conf -> N)*(conf -> N) % conf -> Micro_n_steps == 0 ) valid = 1; 
		else { printf(" The number of subintervals must divide N*N. Try again. \n"); 
		}	
	}
	while(!valid); 
	// Here I need to put a check that Plot_Micro divides N*N 
	// Ask for the number of simulations to be made
	printf(" Choose the number of simulations N_similations to be made \n ");
	scanf( "%d",&conf -> N_simulations); 
	choose_rates( &conf -> annihilation, &conf -> creation, conf -> N); //Beware of the order of the inputs
	choose_initialization_function(conf); 
	}

	void print_config_to_screen(Configuration *conf){	
         printf("You have made the following choices: \n");
	 printf("   N: %d,\n\n",conf ->  N);
   	 printf("   T: %lf \n\n", conf -> T);
   	 printf("   L: %d \n\n", conf -> L);
	 printf("   N_simulations: %d \n\n", conf -> N_simulations);
  	 printf("   Resolution: %d,\n\n", conf -> resolution);
	 printf("   Micro_n_steps: %d \n\n", conf -> Micro_n_steps );
	 printf("   Micro_Creation_Rate: %lf \n\n",conf -> creation );
	 printf("   Micro_Annihilation_Rate: %lf \n\n\n",conf -> annihilation); 
	 printf("   A1_magnetization: %lf \n\n\n", conf -> A1_magnetization); 
	 printf("   B1_magnetization: %lf \n\n\n",conf -> B1_magnetization); 
	 printf("   A1_corr: %lf \n\n\n",conf -> A1_corr); 
	 printf("   A2_corr: %lf \n\n\n",conf -> A2_corr); 
	 printf("   B1_corr: %lf \n\n\n",conf -> B1_corr); 
	 printf("   B2_corr: %lf \n\n\n",conf -> B2_corr); 
	 printf("   A1_time_delayed_corr: %lf \n\n\n",conf -> A1_time_delayed_corr); 
	 printf("   A2_time_delayed_corr: %lf \n\n\n",conf -> A2_time_delayed_corr); 
	 printf("   B1_time_delayed_corr: %lf \n\n\n",conf -> B1_time_delayed_corr); 
	 printf("   B2_time_delayed_corr: %lf \n\n\n",conf -> B2_time_delayed_corr); 
//	 printf("   Lambda_T: %lf \n\n\n",conf -> Lambda_T); 
	}
	
//	printf("From which we obtain: \n" ); 
//	 printf("   l:  %d \n\n\n",conf -> parameters -> l); 
//	 printf("   Total sites :  %d \n\n\n",conf -> parameters -> tot_sites); 
//	 printf("   n_micro_evolution : %d \n\n\n",conf -> parameters -> n_micro_evolution); 
//	 printf("   a1: %d \n\n\n",conf -> corr_par -> a1) ; 
//	 printf("   a2: %d \n\n\n",conf -> corr_par -> a2); 
//	 printf("   b1: %d \n\n\n",conf -> corr_par -> b1); 
//	 printf("   b2: %d \n\n\n",conf -> corr_par -> b2); 
//	}


void read_config_file(Configuration *conf, const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening config file");
        exit(EXIT_FAILURE);
    }

    char line[MAX_LINE], key[64], value[128];
    while (fgets(line, sizeof(line), file)) {
        if (is_comment_or_blank(line)) continue;
        if (sscanf(line, "%63[^=]=%127s", key, value) != 2) continue;

        int found = 0;
        for (int i = 0; config_table[i].key; i++) {
            if (strcmp(config_table[i].key, key) == 0) {
                if (config_table[i].setter(conf, value) != 0) {
                    fprintf(stderr, "Invalid value for key: %s\n", key);
                    exit(EXIT_FAILURE);
                }
                found = 1;
                break;
            }
        }

        if (!found) {
            fprintf(stderr, "Unknown configuration key: %s\n", key);
            exit(EXIT_FAILURE);
        }
    }

    fclose(file);

    choose_rates(&conf->annihilation, &conf->creation, conf->N);
    choose_initialization_function(conf);
 
}



void print_conf_json(const Configuration *conf, const char *dirname){
	 char filepath[512];
    	 snprintf(filepath, sizeof(filepath), "%s/configuration.json", dirname);
	 FILE *conf_json = fopen(filepath, "w");
	  if (!conf_json) {
	 	perror("Failed to write config.json");
	 	exit(EXIT_FAILURE);
	 } 
         fprintf(conf_json, "{\n");
	 fprintf(conf_json, "  \"N\": %d,\n",conf ->  N);
   	 fprintf(conf_json, "  \"T\": %lf, \n", conf -> T);
   	 fprintf(conf_json, "  \"L\": %d, \n", conf -> L);
	 fprintf(conf_json, "  \"N_simulations\": %d, \n", conf -> N_simulations);
  	 fprintf(conf_json, "  \"resolution\": %d,\n", conf -> resolution);
	 fprintf(conf_json, "  \"Micro_n_steps\": %d, \n", conf -> Micro_n_steps );
	 fprintf(conf_json, "  \"Micro_Creation_Rate\": %.5f, \n",conf -> creation );
	 fprintf(conf_json, "  \"Micro_Annihilation_Rate\": %.5f \n",conf -> annihilation);
  	 fprintf(conf_json, "}\n"); 
    	 fclose(conf_json); 
}




void choose_rates(double *annihilation, double *creation, int N ) {	
	char input[10]; 
	int valid = 0; 
	int choice; 
	do {
	printf(" Choose the  creation and the annihilation rates for the microscopic dynamics \n\n " );
	printf(" 1. Creation: 2/(N*N) Annihilation: 2 \n 2. Creation: 2/N  Annihilation: 2N \n 3. Creation: 0, Annihilation: 0 \n" );
	printf(" 4. Creation: 10/ N*N, Annihilation: 2\n 5. Creation: 2/N^2 Annihilation: 20 \n 6. Free choice \n");
	// 1. Is the choice of my master thesis, and is the case in which the generator has a pointwise limit, see the file Glauber.tex 
	// 2. Is the choice we are making 
	// 3. Is the non noisy voter model, and the invariant measure is obtaining by independently placing the spins with +1 or -1  
	// 4. The invariant measure is not the Ising one, but it should be similar	
		clear_stdin();
		if(fgets(input, sizeof(input), stdin)) {
			if(sscanf(input, "%d", &choice )){
				if( choice >= 1 && choice < 67) { 
					valid = 1; 
				} else {
					printf(" Invalid choice, try again \n"); 
				} 
			} else {
				printf(" Invalid input. Please enter a number.\n");	
			}
		}  else{
			printf("Error reading input\n"); 
		}
	}
	while(!valid);
	// Choice 1 are (a slight) modification of the rates chosen by Glauber
	switch( choice ) {

	case 1: {                     
		 *creation = 2.0/((double) N*N); 
		 *annihilation = 2.0;  
		break; 
		}
	// Choice 2 has the same invariant measure, but creation and annihilation are speeded up by a factor $N$
	case 2: {	
		 *creation = 2.0/((double) N);
		 *annihilation = 2.0*((double ) N);         
		break; 
		}
	// Choice 3
	case 3: {
		 *annihilation = 0; 
		 *creation = 0;        
		break; 
	}
	case 4: {
		 *creation = 10.0/((double) N*N); 
		 *annihilation = 1.0;           
		break; 
	}
	case 5: {
		 *creation = 2.0/((double) N*N);  // If the creation rate is 2/N^2 for each microscopic timestep 
		 *annihilation = 20.0;  // Here is 1/N^2, so more interfaces get microscopically created
		break; 
	}
	case 6: {
		printf(" Choose the annihilation rate  \n" ); 
		scanf("%lf", annihilation); 
		while ((getchar()) != '\n');
		printf(" Choose the creation rate  \n"); 
		scanf("%lf", creation); 
		while ((getchar()) != '\n');
		break; 
	}
	}
}



void set_parameters( Configuration  *conf, Parameters_update *parameters ){
	// We save the number of simulations to be made at each macroscopic time step N*N
	// factorized using two numbers, since after each N*N/ we have to plot the microscopic spin histoty
	// We set the number of microscopic sites
	 parameters -> l = (int) (conf -> L)*(conf -> N);  // 
	 parameters -> prob_movement = 1 - exp( - 1.0 /(( (double) conf -> Micro_n_steps)) ); // form the rates we compute the probability of jumping in the interval dT*N*N. The rates should be multiplied by N*N, but, since we are diffusively rescaling the picture. 
	 parameters -> prob_annihilation = 1 - exp( - ((double)conf -> annihilation )/(( (double) conf -> Micro_n_steps))); // Rates of the microscopic process 
	
	 parameters -> prob_creation = 1 - exp( -( (double) conf -> creation) / (( (double) conf -> Micro_n_steps))); 
	
	 parameters -> n_micro_evolution = (int)  2*(conf -> Micro_n_steps )*(conf -> N)*(conf -> N)/(conf -> resolution); // The factor 2 is due to the fact that I am updating only a part of the spins at each Glauber update.  

	////////// MATHEMATICAL DEBUG /////////////////////////////////////////////////////////////////////////////
//	parameters -> prob_annihilation =(conf -> N)* (conf -> N)* (parameters -> prob_creation); 
	// The above position is the right Glauber discrete dynamics. The estimation of the probability made through the rates gives an undervalued annihilation probability. 
	// Here I will write what differences I see

	printf("DEBUG INFORMATION: creation rate =%.6f, Micro_n_steps =%d, creation probability =%.6f\n",
        conf->creation, conf->Micro_n_steps, parameters ->prob_creation);

printf("DEBUG INFORMATION: Annihilation rate =%lf,  Micro_n_steps= %d, annihilation probability =%.6f\n",
       conf -> annihilation,  conf->Micro_n_steps,  parameters -> prob_annihilation);

printf("DEBUG INFORMATION: N^2 * creation probability  =%.6f, prob_annihilation=%.6f\n",
       ((double ) conf -> N) *( (double) conf -> N)* ( parameters -> prob_creation), parameters->prob_annihilation);

// INITIALIZING THE MERSENNE TWISTER

	const gsl_rng_type *G; 
	gsl_rng_env_setup();
	G = gsl_rng_mt19937;
	parameters -> rng = gsl_rng_alloc(G);
}


void choose_initialization_function(Configuration *conf){
	char input[10]; 
	int valid = 0; 
	int choice; 
	do{
		printf(" We now choose the distribution for the initial data for each simulation \n");
		printf(" 1. The initial spin configuration is the Ising one of parameter lambda/N, lambda ro b echosen \n" );
		printf(" 2. We start with two interfaces, one at 0 and the other one at L/2 \n"); 
		printf(" 3. We start without interfaces\n "); 
		printf(" 4. We start with a Bernoulli measure of parameter p \n");
		printf(" 5. We start with a configuration with N  interfaces \n"); 
		fflush(stdout); 
		if(fgets( input, sizeof(input),stdin)){
			if(sscanf(input , "%d", &choice)){
				if( choice < 6 && choice >=1  )  valid = 1;
				else printf("Invalid choice. Please choose a number between 1 and 4"); 
			} else{
				printf("Error. Please type a number \n"); 	
			}
		} else{
			printf("Error in reading the choice\n"); 
		}
	} 
	while(!valid); 
	switch(choice){
	
	case 1: {
			conf -> initialize = initialize_spins_1; 
			double *p = malloc(sizeof(double)); 
			valid = 0; 
			do {	
			printf(" The geometric process has parameter lambda/N. Choose lambda \n " );
			fflush(stdout); 
			clear_stdin(); 
			scanf("%lf" , p );
			if( *p < 0) {
				printf("Invalid choice. Please insert positive number"); 
			} else {
				valid = 1; 
			}
			}
			while(!valid); 	
			(*p) = (*p)/(conf -> N); 
			printf(" The parameter is %f", *p ); 
			conf -> initialization_param = p;
			break; 
		}	
	case 2: {
			conf -> initialize = initialize_spins_2; 
			conf -> initialization_param = NULL; 
			break; 
		}	
	case 3: {
			conf -> initialize = initialize_spins_3; 
			conf -> initialization_param = NULL; 
			break;
		}	
	case 4: {
			conf -> initialize = initialize_spins_4; 
			double *p = malloc(sizeof(double)); 
			valid = 0; 
			clear_stdin();
			do {	
			printf(" Choose the parameeter p of the Bernoullian initial measure\n " ); 
			fflush(stdout);
			scanf("%lf" , p );
			if( *p < 0 || *p> 1) {
				printf("Invalid choice. Please insert a number between 0 and 1"); 
			} else {
				valid = 1; 
			}
			}
			while(!valid); 	
			printf("%f", *p ); 
			conf -> initialization_param = p; 
			break;
		}

	case 5: {
	conf -> initialize = initialize_spins_5;
	conf -> initialization_param = NULL; 
	break;
	}	

//	case 6: {
//			conf -> initialize = initialize_spins_6; 
//			double *p = malloc(sizeof(double)); 
//			valid = 0; 
//			do {	
//			printf(" The geometric process has parameter lambda/N. Choose lambda \n " );
//			fflush(stdout); 
//			clear_stdin(); 
//			scanf("%lf" , p );
//			if( *p < 0) {
//				printf("Invalid choice. Please insert positive number"); 
//			} else {
//				valid = 1; 
//			}
//			}
//			while(!valid); 	
//			(*p) = (*p)/(conf -> N); 
//			printf(" The parameter is %f", *p ); 
//			conf -> initialization_param = p;
//			break; 
//		}	
	}	
}

void initialize_spins_1(int *spins, int l, void *params) {
    spins[0] = (rand() % 2) * 2 - 1; // ±1
    for (int i = 1; i < l; i++) {
        double *p =(double*) params;
        if ((rand() / (double)RAND_MAX) < *p)
            spins[i] = -spins[i - 1];
        else
            spins[i] = spins[i - 1];
    }
}

void initialize_spins_2(int *spins, int  l, void *params) {
    for (int i = 0; i < l; i++) {
    	spins[i] = -1; 
	}
    for (int i = 0; i < l/2; i++ ){
	spins[i] = 1; 
	}
}


void initialize_spins_3(int *spins, int l, void *params) {
    for (int i = 0; i < l; i++) {
    	spins[i] = -1; 
	}
}


void initialize_spins_4(int *spins, int l, void *params) {
    double *p = (double *) params ;
    for(int i = 0; i < l; i++ ){
		if(rand()/ (double) RAND_MAX < *p ){
			spins[i] = 1; 
	} else { 
		spins[i] = -1; 
	}
    }
} 


void initialize_spins_5( int *spins, int l,  void *params ){
	for(int i = 0; i < l ; i ++  ){
		if (i % 2 == 0) spins[i] = 1;
		else{ spins[i] = -1;
		} 
		
	}
}



void initialize_spins_6(int *data, int l,  int k, double p,  int tot_sites) {
    // Preconditions:  l>=2, k<l/2; 
    //  Is the Ising probability conditioned of on having spins +1 from l/2 to l + k 
    data[0] = (rand() % 2) * 2 - 1; //
    for (int i = 1; i < tot_sites; i++) {
		if ( l/2 < i % l <l/2 + k ){
			data[i] = 1; 
		}
		else {
	      	 	if ((rand() / (double)RAND_MAX) < p){
       			     data[i] = -data[i - 1];
      	 		} else{
        		     data[i] = data[i - 1];
	    		}
		}
	}
}

void Glauber_update_parity(int *spins_aux , int *new_spins, Parameters_update *params,int a){ // a can be either 0 or 1, depending on its if I want to update odd or even sites
	int left, right, si; // spin auxiliary variables
	double p;  
	int count = 0;  
	for(int s = 0 ; s < (*params).n_micro_evolution; s++){
	memcpy(spins_aux, new_spins, (params -> l )*(sizeof(int))); 
				if( a == 1 ){
					 p = gsl_rng_uniform((*params).rng);
					 left = spins_aux[(*params).l -1 ];
               				 right = spins_aux[1];
     		  	      		 si = spins_aux[0]; 
					 if (left == right){


                   				 if (left == si){ // we are in the case of an unfavourable spin
                     				   	if (p < (*params).prob_creation){ 
								new_spins[0] = - new_spins[0];
						 	}
						} 
		   				if (left != si ) {
							if ( p <(*params).prob_annihilation) {
								new_spins[0] = - new_spins[0]; 
							}
						}
		   			 } else { 
					        if ( p < (*params).prob_movement ) {
							new_spins[0] = - new_spins[0]; 
						}
					 }
				// Here is the second copy 
				//// I am assuming that l is an even number. Therefore, if I update the spin at 0 , I do not need to update the spin at l-1
				} else{	 p = gsl_rng_uniform((*params).rng);
					 left = spins_aux[(*params).l-2];
               				 right = spins_aux[0];
					si = spins_aux[params -> l -1]; 
		  			 if (left == right){
                   				 if (left == si){ // we are in the case of an unfavourable spin
                     				   	if (p < (*params).prob_creation) {
								new_spins[(*params).l-1] = - new_spins[(*params).l-1];
							 }
						} 
		   				if (left != si ) {
							if ( p <(*params).prob_annihilation){ 
								new_spins[(*params).l-1] = - new_spins[(*params).l-1]; 
							}
						}
		   			 } else { 
					        if ( p < (*params).prob_movement ) {
							new_spins[(*params).l-1] = - new_spins[(*params).l-1]; 
						}
		   			 }	
				}		
			// Here for the remaining sites. We are assuming $l>=2$. 
			for (int i = 1 + a; i < (*params).l - 1 ; i += 2) {
               	                 p = gsl_rng_uniform((*params).rng);
				 left = spins_aux[(i - 1)];
               			 right = spins_aux[(i + 1) ];
       			         si = spins_aux[i]; 
 		  		 if (left == right){
                   			 if (left == si){ // we are in the case of an unfavourable spin
                     			   	if (p < (*params).prob_creation) {
						new_spins[i] = - new_spins[i];
						}
					  } 
		   			if (left != si ) {
						if ( p <(*params).prob_annihilation){ 
							new_spins[i] =  - new_spins[i]; 
						}
					}
		   		 } else { 
				        if ( p < (*params).prob_movement ){ 
						new_spins[i] = - new_spins[i]; 
					}		 
	  		 	}			
			}	
	}
}






//void print_spins_to_file_binary(int *spins,  FILE *f , Parameters *param){
//		fwrite(spins, sizeof(int), param -> l ,f); 
//}








//
//void Glauber_update2(int *spins_aux , int *new_spins, Parameters_update *params){
//	int left, right, si; // spin auxiliary variables
//	double p;   
//	for(int s =0; s < (*params).n_micro_evolution; s++){
//		memcpy(spins_aux, new_spins, (params -> tot_sites)*(sizeof(int))); 
//		for (int i = 0; i < (params -> tot_sites) ; i++) {  
//			p = gsl_rng_uniform((*params).rng);
//			if( i % (params -> l)  == 0 ){
//				 left = spins_aux[i + (params -> l) - 1]; 
//			} else {
//				 left =  spins_aux[i -1];
//			}
//			if (i % (params -> l) == params -> l -1 ){
//				right = spins_aux[i - (params -> l) + 1]; 
//			} else{
//				right = spins_aux[i+1]; 
//			}	
//       			si = spins_aux[i]; 
// 		  	if (left == right){
//                   		 if (left == si){ // we are in the case of an unfavourable spin
//                     			if (p < (*params).prob_creation) {
//						new_spins[i] = - new_spins[i];
//					}
//				 } 
//		   		if (left != si ) {
//					if ( p <(*params).prob_annihilation){ 
//						new_spins[i] =  - new_spins[i]; 
//					}
//				}
//		   		 } else { 
//				        if ( p < (*params).prob_movement ){ 
//						new_spins[i] = - new_spins[i]; 
//					}		 
//	  		 	}			
//			}	
//	}
//}
//


//void Glauber_update(int *spins_aux , int *new_spins, Parameters_update *params){
//	int left, right, si; // spin auxiliary variables
//	double p;  
//	int count = 0;  
//	for(int s =0; s < (*params).n_micro_evolution; s++){
//	memcpy(spins_aux, new_spins, (params -> l )*(sizeof(int))); 
//				 p = gsl_rng_uniform((*params).rng);
//				 left = spins_aux[(*params).l -1 ];
//               			 right = spins_aux[1];
//     		  	      	 si = spins_aux[0]; 
//				 if (left == right){
//
//                   			 if (left == si){ // we are in the case of an unfavourable spin
//                     			   	if (p < (*params).prob_creation){ 
//							new_spins[0] = - new_spins[0];
//					 	}
//					} 
//		   			if (left != si ) {
//						if ( p <(*params).prob_annihilation) {
//							new_spins[0] = - new_spins[0]; 
//						}
//					}
//		   		 } else { 
//				        if ( p < (*params).prob_movement ) {
//						new_spins[0] = - new_spins[0]; 
//					}
//				 }
//				// Here is the second copy 
//				 p = gsl_rng_uniform((*params).rng);
//				 left = spins_aux[(*params).l-2];
//               			 right = spins_aux[0];
//				si = spins_aux[params -> l -1]; 
//		  		 if (left == right){
//                   			 if (left == si){ // we are in the case of an unfavourable spin
//                     			   	if (p < (*params).prob_creation) {
//							new_spins[(*params).l-1] = - new_spins[(*params).l-1];
//						 }
//					} 
//		   			if (left != si ) {
//						if ( p <(*params).prob_annihilation){ 
//							new_spins[(*params).l-1] = - new_spins[(*params).l-1]; 
//						}
//					}
//		   		 } else { 
//				        if ( p < (*params).prob_movement ) {
//						new_spins[(*params).l-1] = - new_spins[(*params).l-1]; 
//					}
//		   		 }			
//			// Here for the remaining sites. We are assuming $l>=2$. 
//			for (int i = 1; i < (*params).l - 1 ; i++) {
//               	                 p = gsl_rng_uniform((*params).rng);
//				 left = spins_aux[(i - 1)];
//               			 right = spins_aux[(i + 1) ];
//       			         si = spins_aux[i]; 
// 		  		 if (left == right){
//                   			 if (left == si){ // we are in the case of an unfavourable spin
//                     			   	if (p < (*params).prob_creation) {
//						new_spins[i] = - new_spins[i];
//						}
//					  } 
//		   			if (left != si ) {
//						if ( p <(*params).prob_annihilation){ 
//							new_spins[i] =  - new_spins[i]; 
//						}
//					}
//		   		 } else { 
//				        if ( p < (*params).prob_movement ){ 
//						new_spins[i] = - new_spins[i]; 
//					}		 
//	  		 	}			
//			}	
//	}
//}
//
//void print_spins_to_file_binary(int *spins,  FILE *f , Parameters_update *param){
//		fwrite(spins, sizeof(int), param -> l ,f); 
//}

void print_spins_to_file_binary2(int *spins,  FILE *f , int length){
	fwrite(spins, sizeof(int), length ,f); 
}


void create_output_directory(const Configuration *conf, char *dirname_out, size_t len) {
    snprintf(dirname_out, len, "results_N%d_T%.1lf_L%d_Annihilation%.2f",
             conf -> N, conf -> T, conf->L, conf->annihilation);

    if (mkdir(dirname_out, 0755) && errno != EEXIST) {
        perror("mkdir failed");
        exit(EXIT_FAILURE);
    }
}


int set_N(Configuration *conf, const char *value) {
      printf("set_T called with value='%s'\n", value);
    int v = atoi(value);
    if (v <= 0) return -1;
    printf(" The value of N: %d", v); 
    conf->N = v;
    printf(" The value of N: %d", conf -> N); 
    return 0;
}

int set_T(Configuration *conf, const char *value) {
    float v = atof(value);
    if (v <= 0.0f) return -1;
    conf->T = v;
    return 0;
}

int set_L(Configuration *conf, const char *value) {
    float v = atof(value);
    if (v <= 0.0f) return -1;
    conf->L = v;
    return 0;
}

int set_N_simulations(Configuration *conf, const char *value) {
    int v = atoi(value);
    if (v <= 0) return -1;
    conf->N_simulations = v;
    return 0;
}

int set_resolution(Configuration *conf, const char *value) {
    int v = atoi(value);
    if (v <= 0) return -1;
    conf->resolution = v;
    return 0;
}

int set_Micro_n_steps(Configuration *conf, const char *value) {
    int v = atoi(value);
    if (v <= 0) return -1;
    conf->Micro_n_steps = v;
    return 0;
}
int set_A1_magnetization(Configuration *conf, const char *value) {
    float v = atof(value);
    conf->A1_magnetization = v;
    return 0;
}

int set_B1_magnetization(Configuration *conf, const char *value) {
    float v = atof(value);
    conf->B1_magnetization = v;
    return 0;
}

int set_A1_corr(Configuration *conf, const char *value) {
    float v = atof(value);
    conf->A1_corr = v;
    return 0;
}

int set_A2_corr(Configuration *conf, const char *value) {
    float v = atof(value);
    conf->A2_corr = v;
    return 0;
}

int set_B1_corr(Configuration *conf, const char *value) {
    float v = atof(value);
    conf->B1_corr = v;
    return 0;
}

int set_B2_corr(Configuration *conf, const char *value) {
    float v = atof(value);
    conf->B2_corr = v;
    return 0;
}

int set_A1_time_delayed_corr(Configuration *conf, const char *value) {
    float v = atof(value);
    conf->A1_time_delayed_corr = v;
    return 0;
}

int set_A2_time_delayed_corr(Configuration *conf, const char *value) {
    float v = atof(value);
    conf->A2_time_delayed_corr = v;
    return 0;
}

int set_B1_time_delayed_corr(Configuration *conf, const char *value) {
    float v = atof(value);
    conf->B1_time_delayed_corr = v;
    return 0;
}

int set_B2_time_delayed_corr(Configuration *conf, const char *value) {
    float v = atof(value);
    conf->B2_time_delayed_corr = v;
    return 0;
}

//int set_Lambda_T(Configuration *conf, const char *value) {
//    float v = atof(value);
//    conf->Lambda_T = v;
//    return 0;
//}


void validate_config(const Configuration *conf) {
    if ((conf->N * conf->N * conf->Micro_n_steps) % conf->resolution != 0) {
        fprintf(stderr, "Micro_n_steps must satisfy N*N*Micro_n_steps divisible by resolution\n");
        exit(EXIT_FAILURE);
    }

    float coords[] = {
        conf->A1_magnetization, conf->B1_magnetization,
        conf->A1_corr, conf->A2_corr, conf->B1_corr, conf->B2_corr,
        conf->A1_time_delayed_corr, conf->A2_time_delayed_corr,
        conf->B1_time_delayed_corr, conf->B2_time_delayed_corr
    };

    for (int i = 0; i < sizeof(coords)/sizeof(float); i++) {
        if (coords[i] <= 0 || coords[i] >= conf->L || !is_scaled_integer(coords[i], conf->N)) {
            fprintf(stderr, "Statistic parameter %.2f is invalid: must be >0, <L, and N*value an integer\n", coords[i]);
            exit(EXIT_FAILURE);
        }
    }
}
int is_comment_or_blank(const char *line) {
    // Skip leading whitespace
    while (isspace(*line)) line++;

    return (*line == '#' || *line == '/' || *line == '\0' || *line == '\n');
}


