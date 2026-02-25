# Lean 4 Setup and Usage Guide

## What is Lean?

Lean 4 is a functional programming language and interactive theorem prover. Lean files (`.lean`) are:
- **Text files** written in the Lean language
- **Not generated** - you write them directly (like Python or C++ files)
- **Checked/verified** by the Lean compiler and proof assistant

## How Lean Files Work

### 1. Lean Files are Source Code

Lean files are just text files with `.lean` extension. You write them like any other code:

```lean
-- This is a comment
def my_function (x : ℕ) : ℕ := x + 1

theorem my_theorem (x : ℕ) : my_function x = x + 1 := by rfl
```

### 2. Lean Project Structure

A Lean project typically has this structure:

```
my_project/
├── lean-toolchain          # Specifies Lean version
├── lakefile.lean          # Build configuration (like Makefile)
└── MyProject/
    └── Basic.lean         # Your Lean source files
```

## Setting Up Lean 4

### Option 1: Using VS Code (Recommended)

1. **Install VS Code**: https://code.visualstudio.com/

2. **Install Lean 4 Extension**:
   - Open VS Code
   - Go to Extensions (Ctrl+Shift+X)
   - Search for "Lean 4"
   - Install the official "Lean 4" extension by leanprover

3. **Install Lean 4**:
   - The extension will prompt you to install Lean 4
   - Or install manually: https://leanprover-community.github.io/get_started.html

4. **Open Your Project**:
   - Open the folder containing your `.lean` files in VS Code
   - The extension will automatically detect and set up the project

### Option 2: Command Line Installation

1. **Install elan** (Lean version manager):
   ```bash
   curl https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -sSf | sh
   ```

2. **Install Lean 4**:
   ```bash
   elan toolchain install stable
   elan default stable
   ```

3. **Verify installation**:
   ```bash
   lean --version
   ```

## Creating a Lean Project

### For Your Biclique Project

1. **Create project structure**:
   ```bash
   cd /data1/yizhangh/ldp-pq
   mkdir -p BicliqueProject
   cd BicliqueProject
   ```

2. **Initialize Lean project**:
   ```bash
   lake init BicliqueProject
   ```

3. **Create lean-toolchain file** (if not created):
   ```bash
   echo "leanprover/lean4:stable" > lean-toolchain
   ```

4. **Copy your Lean file**:
   ```bash
   cp ../ButterflyCounting_2_2.lean BicliqueProject/BicliqueProject/ButterflyCounting.lean
   ```

## Running/Checking Lean Files

### Method 1: VS Code (Interactive)

1. **Open the file** in VS Code
2. **Lean server starts automatically**
3. **See real-time feedback**:
   - Green checkmark ✓ = proof is correct
   - Red error × = error in code
   - Yellow warning ⚠ = incomplete proof (has `sorry`)

4. **Interact with proofs**:
   - Hover over code to see types
   - Click on `sorry` to see what needs to be proved
   - Use "Go to Definition" (F12) to navigate

### Method 2: Command Line

1. **Check a single file**:
   ```bash
   lean --check ButterflyCounting_2_2.lean
   ```

2. **Build the project**:
   ```bash
   lake build
   ```

3. **Run in interactive mode**:
   ```bash
   lean --server
   ```

### Method 3: Using Lake (Build System)

```bash
# Build the project
lake build

# Check all files
lake env lean --check BicliqueProject/ButterflyCounting.lean

# Run tests (if you have them)
lake test
```

## Understanding Lean Output

### Successful Check
```
ButterflyCounting_2_2.lean:314:0: information: declaration uses 'sorry'
```
This means the file compiles but has incomplete proofs.

### Error Example
```
ButterflyCounting_2_2.lean:45:15: error: unknown identifier 'some_function'
```
This means there's a syntax or type error.

### Successful Proof
When a proof is complete (no `sorry`), you'll see:
- No errors
- Green checkmark in VS Code
- File compiles successfully

## Working with Your Butterfly File

### Current Status

Your `ButterflyCounting_2_2.lean` file:
- ✅ **Compiles** (no syntax errors)
- ⚠️ **Has `sorry`** (incomplete proofs)
- 📝 **Ready for proof development**

### Next Steps to Complete Proofs

1. **Open in VS Code**:
   ```bash
   code /data1/yizhangh/ldp-pq/ButterflyCounting_2_2.lean
   ```

2. **See what needs proving**:
   - Look for lines with `sorry`
   - Hover to see the goal

3. **Start proving**:
   - Replace `sorry` with actual proof tactics
   - Lean will check in real-time

### Example: Completing a Proof

**Before** (with `sorry`):
```lean
theorem example_theorem (x : ℕ) : x + 0 = x := by
  sorry
```

**After** (complete proof):
```lean
theorem example_theorem (x : ℕ) : x + 0 = x := by
  rfl  -- "reflexivity" - proves equality by definition
```

## Lean Project Setup for Your Files

### Quick Setup Script

Create `setup_lean_project.sh`:

```bash
#!/bin/bash
# Setup Lean project for biclique verification

cd /data1/yizhangh/ldp-pq

# Create project directory
mkdir -p BicliqueLean
cd BicliqueLean

# Initialize Lean project
lake init BicliqueLean

# Create lean-toolchain
echo "leanprover/lean4:stable" > lean-toolchain

# Copy Lean files
mkdir -p BicliqueLean
cp ../ButterflyCounting_2_2.lean BicliqueLean/ButterflyCounting.lean

# Update lakefile.lean to include Mathlib
# (You'll need to add Mathlib dependencies)

echo "Project created! Run 'lake build' to build."
```

### Adding Mathlib Dependencies

Edit `lakefile.lean`:

```lean
import Lake
open Lake DSL

package «BicliqueLean» where
  -- Add your package configuration here

@[default_target]
lean_lib «BicliqueLean» where
  -- Add your library configuration here

-- Add Mathlib dependency
require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "master"
```

Then run:
```bash
lake update
lake build
```

## Common Commands

### Development Workflow

```bash
# 1. Check if file compiles
lean --check ButterflyCounting_2_2.lean

# 2. Build project
lake build

# 3. Update dependencies
lake update

# 4. Clean build artifacts
lake clean

# 5. Run in interactive mode (for debugging)
lean --server
```

### VS Code Shortcuts

- **Ctrl+Shift+P** → "Lean 4: Restart Server" (if things get stuck)
- **Hover** over code → See types and documentation
- **F12** → Go to definition
- **Ctrl+Click** → Jump to definition
- **Alt+Enter** → Apply suggested fix

## Troubleshooting

### Issue: "Lean server not responding"

**Solution**:
1. Restart Lean server: Ctrl+Shift+P → "Lean 4: Restart Server"
2. Check Lean is installed: `lean --version`
3. Rebuild project: `lake clean && lake build`

### Issue: "Unknown identifier"

**Solution**:
1. Check imports at top of file
2. Make sure Mathlib is installed: `lake update`
3. Check spelling of identifiers

### Issue: "Cannot find declaration"

**Solution**:
1. The declaration might be in a different file - check imports
2. Might need to add to `lakefile.lean`
3. Run `lake build` to compile dependencies

### Issue: File won't compile

**Solution**:
1. Check syntax (missing colons, parentheses, etc.)
2. Check all imports are available
3. Look at error messages - they're usually helpful
4. Try commenting out problematic sections

## Example: Checking Your File

```bash
# Navigate to your project
cd /data1/yizhangh/ldp-pq

# Check the file directly (if standalone)
lean --check ButterflyCounting_2_2.lean

# Or set up as proper project
mkdir -p BicliqueLean
cd BicliqueLean
lake init BicliqueLean
cp ../ButterflyCounting_2_2.lean BicliqueLean/ButterflyCounting.lean
lake build
```

## Resources

- **Lean 4 Manual**: https://leanprover.github.io/lean4/doc/
- **Mathlib Documentation**: https://leanprover-community.github.io/
- **Lean 4 Zulip Chat**: https://leanprover.zulipchat.com/ (for questions)
- **Theorem Proving in Lean 4**: https://leanprover.github.io/theorem_proving_in_lean4/ (tutorial)

## Summary

1. **Lean files are text files** - you write them directly
2. **Install Lean 4** via VS Code extension or command line
3. **Open in VS Code** for best experience (real-time checking)
4. **Use `lean --check`** for command-line verification
5. **Use `lake build`** for project builds
6. **Replace `sorry`** with actual proofs to complete verification

Your `ButterflyCounting_2_2.lean` file is ready to use - just open it in VS Code with the Lean extension installed!
