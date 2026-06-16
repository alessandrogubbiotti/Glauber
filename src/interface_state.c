#include "interface_state.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_randist.h>

#define INF 1e300

/* =========================================================================
 * Indexed binary min-heap over the k domain events + 1 creation event.
 * Replaces the O(k) next-event linear scan with O(1) peek / O(log k) update.
 * Node ids: 0..k-1 = domains, CREATE_ID (-1) = the global creation event.
 * The heap is kept valid at every step: each single tau change is followed
 * immediately by heap_update() (correct by the standard decrease/increase-key
 * argument), while structural events (creation, annihilation, wraparound
 * movement) rebuild it in O(k) via heap_build().
 * ========================================================================= */

#define CREATE_ID (-1)

static inline double heap_key(const InterfaceState *s, int id)
{
    return (id == CREATE_ID) ? s->create_tau : s->domain_tau[id];
}

/* Place node `id` at heap slot `pos`, updating its recorded position. */
static inline void heap_place(InterfaceState *s, int pos, int id)
{
    s->heap[pos] = id;
    if (id == CREATE_ID) s->heap_create_pos = pos;
    else                 s->heap_pos[id]    = pos;
}

static void heap_sift_up(InterfaceState *s, int pos)
{
    int    id  = s->heap[pos];
    double key = heap_key(s, id);
    while (pos > 0) {
        int parent = (pos - 1) / 2;
        int pid    = s->heap[parent];
        if (heap_key(s, pid) <= key) break;
        heap_place(s, pos, pid);
        pos = parent;
    }
    heap_place(s, pos, id);
}

static void heap_sift_down(InterfaceState *s, int pos)
{
    int    n   = s->heap_n;
    int    id  = s->heap[pos];
    double key = heap_key(s, id);
    while (1) {
        int    sm  = pos;
        double smk = key;
        int    lft = 2 * pos + 1, rgt = 2 * pos + 2;
        if (lft < n) { double k = heap_key(s, s->heap[lft]); if (k < smk) { sm = lft; smk = k; } }
        if (rgt < n) { double k = heap_key(s, s->heap[rgt]); if (k < smk) { sm = rgt; smk = k; } }
        if (sm == pos) break;
        heap_place(s, pos, s->heap[sm]);
        pos = sm;
    }
    heap_place(s, pos, id);
}

/* One node's key changed; restore the heap (the rest is assumed valid). */
static inline void heap_update(InterfaceState *s, int id)
{
    int pos = (id == CREATE_ID) ? s->heap_create_pos : s->heap_pos[id];
    heap_sift_up(s, pos);
    pos = (id == CREATE_ID) ? s->heap_create_pos : s->heap_pos[id];
    heap_sift_down(s, pos);
}

/* Rebuild the whole heap from the current domain_tau[] and create_tau (O(k)).
 * Used after structural events that add/remove domains or reindex the array. */
static void heap_build(InterfaceState *s)
{
    s->heap_n = s->k + 1;
    for (int j = 0; j < s->k; j++) heap_place(s, j, j);
    heap_place(s, s->k, CREATE_ID);
    for (int pos = s->heap_n / 2 - 1; pos >= 0; pos--)
        heap_sift_down(s, pos);
}

/* =========================================================================
 * Internal helpers — not exposed in the header.
 * ========================================================================= */

static inline int domain_width(const InterfaceState *s, int j)
{
    int next = (j + 1) % s->k;
    return (s->ifaces[next].pos - s->ifaces[j].pos + s->l) % s->l;
}

/* Rate of the domain-level event for domain j. */
static inline double domain_rate(int w, double ann_rate)
{
    if (w == 1) return ann_rate;
    if (w >= 2) return 2.0;
    return 0.0; /* w==0: impossible during a valid run */
}

/* Bulk spins contributed by a domain of width w. */
static inline int bulk_of(int w)
{
    return (w >= 3) ? w - 2 : 0;
}

/* Record that the spin at site i flipped (no-op unless tracking enabled). */
static inline void mark_flip(InterfaceState *s, int i)
{
    if (s->flipped) s->flipped[i] = 1;
}

/* Sample a fresh absolute event time for domain j from the current micro_time. */
static double sample_domain_tau(InterfaceState *s, int j)
{
    int w    = domain_width(s, j);
    double r = domain_rate(w, s->ann_rate);
    if (r <= 0.0) return INF;
    return s->micro_time + gsl_ran_exponential(s->rng, 1.0 / r);
}

/* Sample a fresh absolute creation time. */
static double sample_create_tau(InterfaceState *s)
{
    double r = (double)s->n_bulk * s->create_rate;
    if (r <= 0.0) return INF;
    return s->micro_time + gsl_ran_exponential(s->rng, 1.0 / r);
}

/* Rescale an existing absolute event time when the underlying rate changes.
 * Uses the NRM memoryless property:
 *   remaining_new = remaining_old * old_rate / new_rate
 * If old_rate == 0 (event didn't exist), sample fresh.
 * If new_rate == 0, the event no longer exists: return INF. */
static double rescale_tau(double old_tau, double old_rate,
                          double new_rate, double t_now)
{
    if (new_rate <= 0.0) return INF;
    if (old_rate <= 0.0 || old_tau >= INF)
        /* Event was absent; we'd need a fresh sample, but caller handles this. */
        return INF; /* caller should sample fresh */
    double remaining_new = (old_tau - t_now) * old_rate / new_rate;
    return t_now + remaining_new;
}

/* Adjust n_bulk and return the delta for update_create_tau purposes. */
static int update_n_bulk(InterfaceState *s, int old_w, int new_w)
{
    int delta = bulk_of(new_w) - bulk_of(old_w);
    s->n_bulk += delta;
    return delta;
}

/* After n_bulk has changed, rescale or resample create_tau.
 * old_create_rate: the rate before the change (= old_n_bulk * create_rate). */
static void update_create_tau(InterfaceState *s, double old_create_rate)
{
    double new_create_rate = (double)s->n_bulk * s->create_rate;
    if (old_create_rate <= 0.0 || s->create_tau >= INF)
        s->create_tau = sample_create_tau(s);
    else
        s->create_tau = rescale_tau(s->create_tau, old_create_rate,
                                    new_create_rate, s->micro_time);
}

/* =========================================================================
 * Insert / remove interfaces (maintain sorted order, shift array).
 * Caller updates domain_tau, n_bulk, etc. afterwards.
 * ========================================================================= */

/* Insert a new interface at bond position pos with given left_color.
 * Returns the index where it was inserted. */
static int insert_interface(InterfaceState *s, int pos, int left_color)
{
    /* Grow array if needed. */
    if (s->k == s->k_alloc) {
        s->k_alloc = s->k_alloc ? s->k_alloc * 2 : 4;
        s->ifaces     = realloc(s->ifaces,     s->k_alloc * sizeof(Interface));
        s->domain_tau = realloc(s->domain_tau, s->k_alloc * sizeof(double));
        s->heap       = realloc(s->heap,       (s->k_alloc + 1) * sizeof(int));
        s->heap_pos   = realloc(s->heap_pos,   s->k_alloc * sizeof(int));
        if (!s->ifaces || !s->domain_tau || !s->heap || !s->heap_pos) {
            perror("insert_interface: realloc failed");
            exit(EXIT_FAILURE);
        }
    }

    /* Find insertion index (sorted by pos). */
    int idx = 0;
    while (idx < s->k && s->ifaces[idx].pos < pos) idx++;

    /* Shift right. */
    memmove(&s->ifaces[idx + 1],     &s->ifaces[idx],
            (s->k - idx) * sizeof(Interface));
    memmove(&s->domain_tau[idx + 1], &s->domain_tau[idx],
            (s->k - idx) * sizeof(double));

    s->ifaces[idx].pos        = pos;
    s->ifaces[idx].left_color = left_color;
    s->domain_tau[idx]        = INF; /* caller must set this */
    s->k++;
    return idx;
}

/* Remove the interface at index j. */
static void remove_interface(InterfaceState *s, int j)
{
    assert(j >= 0 && j < s->k);
    memmove(&s->ifaces[j],     &s->ifaces[j + 1],
            (s->k - j - 1) * sizeof(Interface));
    memmove(&s->domain_tau[j], &s->domain_tau[j + 1],
            (s->k - j - 1) * sizeof(double));
    s->k--;
}

/* =========================================================================
 * Invariant checker (enabled in debug builds).
 * ========================================================================= */

#ifndef NDEBUG
static void assert_invariants(const InterfaceState *s)
{
    assert(s->k % 2 == 0); /* torus: must be even */

    /* k==0: all sites are bulk. */
    if (s->k == 0) {
        assert(s->n_bulk == s->l);
        return;
    }

    int nb = 0;
    for (int j = 0; j < s->k; j++) {
        int next = (j + 1) % s->k;
        assert(s->ifaces[next].left_color == -s->ifaces[j].left_color);
        nb += bulk_of(domain_width(s, j));
    }
    assert(s->n_bulk == nb);

    /* Heap consistency: positions round-trip and the min-heap order holds. */
    if (s->heap) {
        assert(s->heap_n == s->k + 1);
        for (int p = 0; p < s->heap_n; p++) {
            int id = s->heap[p];
            int recorded = (id == CREATE_ID) ? s->heap_create_pos : s->heap_pos[id];
            assert(recorded == p);
            int lft = 2 * p + 1, rgt = 2 * p + 2;
            if (lft < s->heap_n) assert(heap_key(s, id) <= heap_key(s, s->heap[lft]));
            if (rgt < s->heap_n) assert(heap_key(s, id) <= heap_key(s, s->heap[rgt]));
        }
    }
}
#else
static inline void assert_invariants(const InterfaceState *s) { (void)s; }
#endif

/* =========================================================================
 * Public: istate_alloc
 * ========================================================================= */

InterfaceState *istate_alloc(const int *spins, int l, int N,
                             double ann_rate, double create_rate,
                             unsigned long seed)
{
    InterfaceState *s = malloc(sizeof(InterfaceState));
    if (!s) { perror("istate_alloc"); exit(EXIT_FAILURE); }

    s->l           = l;
    s->N           = N;
    s->ann_rate    = ann_rate;
    s->create_rate = create_rate;
    s->micro_time  = 0.0;
    s->k           = 0;
    s->k_alloc     = 0;
    s->n_bulk      = 0;
    s->n_events    = 0;
    s->ifaces      = NULL;
    s->domain_tau  = NULL;
    s->create_tau  = INF;
    s->heap        = NULL;
    s->heap_pos    = NULL;
    s->heap_create_pos = 0;
    s->heap_n      = 0;
    s->flipped     = NULL;   /* tracking off unless istate_track_flips() */
    s->uniform_color = spins[0]; /* exact when k==0; harmless default else */

    s->rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(s->rng, seed);

    /* Scan spins left to right; detect interfaces at bonds. */
    for (int i = 0; i < l; i++) {
        if (spins[i] != spins[(i + 1) % l]) {
            int idx = insert_interface(s, i, spins[i]);
            s->domain_tau[idx] = INF; /* will be set below */
        }
    }

    if (s->k % 2 != 0) {
        fprintf(stderr, "istate_alloc: odd number of interfaces (%d) on torus\n",
                s->k);
        exit(EXIT_FAILURE);
    }

    /* Compute n_bulk and initial domain_tau. */
    s->n_bulk = 0;
    for (int j = 0; j < s->k; j++) {
        s->n_bulk += bulk_of(domain_width(s, j));
        s->domain_tau[j] = sample_domain_tau(s, j);
    }

    /* Handle k==0: entire torus is one domain. */
    if (s->k == 0) {
        /* No domain events; only global creation (all l sites are bulk). */
        s->n_bulk = l; /* all sites are bulk when there are no interfaces */
    }

    s->create_tau = sample_create_tau(s);

    /* If k == 0 no interface was inserted, so insert_interface never allocated
     * the heap arrays; reserve room for the lone creation node. */
    if (!s->heap) s->heap = malloc(1 * sizeof(int));
    heap_build(s);

    assert_invariants(s);
    return s;
}

/* =========================================================================
 * Public: istate_free
 * ========================================================================= */

void istate_free(InterfaceState *s)
{
    gsl_rng_free(s->rng);
    free(s->ifaces);
    free(s->domain_tau);
    free(s->heap);
    free(s->heap_pos);
    free(s->flipped);
    free(s);
}

/* =========================================================================
 * Public: per-site flip tracking (persistence observable)
 * ========================================================================= */

void istate_track_flips(InterfaceState *s)
{
    free(s->flipped);
    s->flipped = calloc((size_t)s->l, sizeof(unsigned char));
    if (!s->flipped) { perror("istate_track_flips"); exit(EXIT_FAILURE); }
}

double istate_never_flipped_fraction(const InterfaceState *s)
{
    assert(s->flipped);
    long n = 0;
    for (int i = 0; i < s->l; i++)
        n += s->flipped[i] ? 0 : 1;
    return (double)n / (double)s->l;
}

/* =========================================================================
 * Public: istate_to_spins
 * ========================================================================= */

void istate_to_spins(const InterfaceState *s, int *spins_out)
{
    if (s->k == 0) {
        /* Uniform state: use the color recorded when the last interface
         * pair annihilated (or at initialisation).  Forcing +1 here would
         * silently flip the sign of every cross-time spin product in
         * realisations that pass through k == 0. */
        for (int i = 0; i < s->l; i++) spins_out[i] = s->uniform_color;
        return;
    }

    /* Fill domains left to right starting from interface 0.
     * Domain j occupies sites ifaces[j].pos+1 .. ifaces[(j+1)%k].pos (wrapping). */
    for (int j = 0; j < s->k; j++) {
        int color = s->ifaces[j].left_color; /* color to the RIGHT of interface j is -left_color */
        int right_color = -color;            /* domain j is to the RIGHT of interface j */

        int start = (s->ifaces[j].pos + 1) % s->l;
        int end   = s->ifaces[(j + 1) % s->k].pos; /* bond position of next interface */

        /* end is the LAST site of domain j (the site to the left of the next interface bond).
         * Sites: start, start+1, ..., end (wrapping on torus). */
        int site = start;
        while (1) {
            spins_out[site] = right_color;
            if (site == end) break;
            site = (site + 1) % s->l;
        }
    }
}

/* =========================================================================
 * Event execution helpers (static)
 * ========================================================================= */

/* Reinsert interface at index j after its pos changed and may have violated
 * sorted order.  Handles the wraparound case (bond l-1 → 0 or 0 → l-1).
 * Returns the new index of the interface. */
static int reinsert_sorted(InterfaceState *s, int j)
{
    Interface iface = s->ifaces[j];
    /* Remove (shifts subsequent elements left, decrements k). */
    memmove(&s->ifaces[j],     &s->ifaces[j + 1],
            (s->k - j - 1) * sizeof(Interface));
    memmove(&s->domain_tau[j], &s->domain_tau[j + 1],
            (s->k - j - 1) * sizeof(double));
    s->k--;
    /* Insert at the correct sorted position (increments k, sets tau=INF). */
    return insert_interface(s, iface.pos, iface.left_color);
}

/* Recompute n_bulk and all domain_tau from scratch (O(k)).
 * Used after wraparound movements where the full state changes. */
static void recompute_all(InterfaceState *s, double old_create_rate)
{
    s->n_bulk = 0;
    for (int jj = 0; jj < s->k; jj++) {
        s->n_bulk       += bulk_of(domain_width(s, jj));
        s->domain_tau[jj] = sample_domain_tau(s, jj);
    }
    update_create_tau(s, old_create_rate);
}

/* execute_movement: domain j fires with w_j >= 2.
 * Direction 0: left interface (j) moves right — domain j shrinks from the left.
 * Direction 1: right interface (j+1)%k moves left — domain j shrinks from the right. */
static void execute_movement(InterfaceState *s, int j)
{
    int dir = (gsl_rng_uniform(s->rng) < 0.5) ? 0 : 1;
    double old_create_rate = (double)s->n_bulk * s->create_rate;

    if (dir == 0) {
        /* Interface j moves RIGHT. */
        int jleft   = (j - 1 + s->k) % s->k;
        int old_wj  = domain_width(s, j);
        int old_wl  = domain_width(s, jleft);
        int old_pos = s->ifaces[j].pos;
        int new_pos = (old_pos + 1) % s->l;
        s->ifaces[j].pos = new_pos;
        mark_flip(s, new_pos);   /* site between bonds old_pos and new_pos */

        if (new_pos < old_pos) {
            /* Wraparound: bond l-1 → 0. Re-sort and recompute everything. */
            reinsert_sorted(s, j);
            recompute_all(s, old_create_rate);
            heap_build(s);
        } else {
            int new_wj = old_wj - 1;
            int new_wl = old_wl + 1;
            update_n_bulk(s, old_wj, new_wj);
            update_n_bulk(s, old_wl, new_wl);
            update_create_tau(s, old_create_rate);
            heap_update(s, CREATE_ID);
            s->domain_tau[j] = sample_domain_tau(s, j);
            heap_update(s, j);
            if (domain_rate(old_wl, s->ann_rate) != domain_rate(new_wl, s->ann_rate)) {
                s->domain_tau[jleft] = sample_domain_tau(s, jleft);
                heap_update(s, jleft);
            }
        }

    } else {
        /* Interface jright = (j+1)%k moves LEFT. */
        int jright  = (j + 1) % s->k;
        int old_wj  = domain_width(s, j);
        int old_wr  = domain_width(s, jright);
        int old_pos = s->ifaces[jright].pos;
        int new_pos = (old_pos - 1 + s->l) % s->l;
        s->ifaces[jright].pos = new_pos;
        mark_flip(s, old_pos);   /* site between bonds new_pos and old_pos */

        if (new_pos > old_pos) {
            /* Wraparound: bond 0 → l-1. Re-sort and recompute everything. */
            reinsert_sorted(s, jright);
            recompute_all(s, old_create_rate);
            heap_build(s);
        } else {
            int new_wj = old_wj - 1;
            int new_wr = old_wr + 1;
            update_n_bulk(s, old_wj, new_wj);
            update_n_bulk(s, old_wr, new_wr);
            update_create_tau(s, old_create_rate);
            heap_update(s, CREATE_ID);
            s->domain_tau[j] = sample_domain_tau(s, j);
            heap_update(s, j);
            if (domain_rate(old_wr, s->ann_rate) != domain_rate(new_wr, s->ann_rate)) {
                s->domain_tau[jright] = sample_domain_tau(s, jright);
                heap_update(s, jright);
            }
        }
    }
}

/* execute_annihilation: domain j fires with w_j == 1.
 * Interfaces j and (j+1)%k are removed; surrounding domains merge. */
static void execute_annihilation(InterfaceState *s, int j)
{
    assert(s->k >= 2);
    assert(domain_width(s, j) == 1);

    /* The lone dissenting spin of the width-1 domain flips. */
    mark_flip(s, (s->ifaces[j].pos + 1) % s->l);

    double old_create_rate = (double)s->n_bulk * s->create_rate;

    /* k==2: both interfaces disappear entirely — the entire torus becomes bulk.
     * The surviving color is that of the surrounding domain, i.e. the color
     * to the LEFT of interface j (the width-1 domain itself has -left_color). */
    if (s->k == 2) {
        s->uniform_color = s->ifaces[j].left_color;
        s->n_bulk = s->l;
        remove_interface(s, 1);
        remove_interface(s, 0);
        /* k is now 0; no domain events remain. */
        update_create_tau(s, old_create_rate);
        heap_build(s);
        return;
    }

    /* General case: k >= 4.
     * jleft is the domain to the LEFT  of the width-1 domain j.
     * jnext is the domain to the RIGHT of the width-1 domain j.
     * jleft != jnext because k >= 4. */
    int jnext = (j + 1) % s->k;
    int jleft = (j - 1 + s->k) % s->k;

    int w_left  = domain_width(s, jleft);
    int w_right = domain_width(s, jnext);
    int w_new   = w_left + 1 + w_right;

    update_n_bulk(s, w_left,  0);
    update_n_bulk(s, 1,       0); /* the width-1 domain */
    update_n_bulk(s, w_right, 0);
    update_n_bulk(s, 0,    w_new);

    /* Remove the two interfaces: higher index first. */
    int hi = (j > jnext) ? j : jnext;
    int lo = (j < jnext) ? j : jnext;
    remove_interface(s, hi);
    remove_interface(s, lo);
    /* s->k is now k-2. */

    /* After the two removals, the merged domain is to the right of the interface
     * that was at index lo-1 (before removal).  In the new array:
     *   lo > 0  → merged domain index = lo - 1
     *   lo == 0 → merged domain is the wrap-around domain = new s->k - 1  */
    int j_merged = (lo == 0) ? (s->k - 1) : (lo - 1);
    s->domain_tau[j_merged] = sample_domain_tau(s, j_merged);

    update_create_tau(s, old_create_rate);
    heap_build(s);
}

/* execute_creation: global creation event fires.
 * Pick a bulk spin uniformly, insert two new interfaces.
 * After insertion all domain taus are recomputed from scratch — this is O(k)
 * but creation events are rare (rate ~ creation per bulk spin), so it is fine. */
static void execute_creation(InterfaceState *s)
{
    assert(s->n_bulk > 0);

    /* Select bulk spin index u in [0, n_bulk). */
    int u = (int)(gsl_rng_uniform(s->rng) * s->n_bulk);
    if (u >= s->n_bulk) u = s->n_bulk - 1; /* clamp floating-point edge */

    int p_left, p_right, domain_color;

    if (s->k == 0) {
        /* All l sites are bulk; pick any site. */
        p_right      = u % s->l;
        p_left       = (p_right - 1 + s->l) % s->l;
        domain_color = s->uniform_color; /* color of the uniform background */
    } else {
        /* Find which domain contains the u-th bulk spin. */
        int chosen_j = -1, local_idx = -1;
        for (int j = 0; j < s->k; j++) {
            int b = bulk_of(domain_width(s, j));
            if (u < b) { chosen_j = j; local_idx = u; break; }
            u -= b;
        }
        assert(chosen_j >= 0);

        /* The bulk spin is at position:
         *   (pos[chosen_j] + 1 + local_idx + 1) % l
         * = (pos[chosen_j] + local_idx + 2) % l
         * New interfaces at bonds p_left = pos+local_idx+1 and p_right = pos+local_idx+2. */
        p_left       = (s->ifaces[chosen_j].pos + local_idx + 1) % s->l;
        p_right      = (s->ifaces[chosen_j].pos + local_idx + 2) % s->l;
        domain_color = -(s->ifaces[chosen_j].left_color); /* color of domain chosen_j */
    }

    /* The nucleated spin (site between bonds p_left and p_right) flips. */
    mark_flip(s, p_right);

    /* Insert the two new interfaces (insert_interface keeps sorted order). */
    insert_interface(s, p_left,  domain_color);
    insert_interface(s, p_right, -domain_color);
    /* k is now k+2. */

    /* Recompute n_bulk and all domain taus from scratch. */
    s->n_bulk = 0;
    for (int j = 0; j < s->k; j++) {
        s->n_bulk       += bulk_of(domain_width(s, j));
        s->domain_tau[j] = sample_domain_tau(s, j);
    }
    s->create_tau = sample_create_tau(s);
    heap_build(s);
}

/* =========================================================================
 * Public: istate_advance_to
 * ========================================================================= */

void istate_advance_to(InterfaceState *s, double target_macro_time)
{
    double target_micro = target_macro_time * (double)s->N * (double)s->N;

    while (s->micro_time < target_micro) {

        /* Next event = heap minimum (O(1) peek).  best_j == -1 → creation. */
        int    best_id  = s->heap[0];
        double best_tau = heap_key(s, best_id);
        int    best_j   = (best_id == CREATE_ID) ? -1 : best_id;

        /* If no event exists (all INF), time cannot advance — absorbing state. */
        if (best_tau >= INF) {
            s->micro_time = target_micro;
            return;
        }

        /* Event falls after the snapshot — do not execute; clamp time. */
        if (best_tau > target_micro) {
            s->micro_time = target_micro;
            return;
        }

        /* Advance time to this event. */
        s->micro_time = best_tau;
        s->n_events++;

        if (best_j < 0) {
            /* Global creation event. */
            s->create_tau = INF; /* consumed; execute_creation rebuilds the heap */
            execute_creation(s);
        } else {
            int w = domain_width(s, best_j);
            s->domain_tau[best_j] = INF; /* consumed */
            heap_update(s, best_j);      /* sink it; keeps the heap valid */
            if (w == 1)
                execute_annihilation(s, best_j);
            else
                execute_movement(s, best_j);
        }

        assert_invariants(s);
    }
}
