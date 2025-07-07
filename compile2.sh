
#!/bin/bash

# Set source file and output name
SRC="Glauber_tidy.c"
OUT="Glauber"

# Optional: Set custom GSL path (e.g. Homebrew on macOS ARM64)
# Uncomment and adjust if needed
# export C_INCLUDE_PATH=/opt/homebrew/include
# export LIBRARY_PATH=/opt/homebrew/lib

# Compile
echo "Compiling $SRC -> $OUT"
#clang "$SRC" -o "$OUT" -lgsl -lgslcblas -lm
clang -g "$SRC" -o "$OUT" -lgsl -lgslcblas -lm
# Check result
if [ $? -eq 0 ]; then
    echo "✅ Build succeeded: ./$OUT"
else
    echo "❌ Build failed"
fi

