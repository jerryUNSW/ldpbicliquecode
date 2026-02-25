# Lean 4 Quick Reference

## TL;DR: How Lean Files Work

1. **Lean files are text files** (`.lean` extension) - you write them like code
2. **They're NOT generated** - you create them manually or with an editor
3. **Lean checks them** - either in VS Code (real-time) or via command line

## Quick Start

### Install Lean 4

**Option 1: VS Code (Easiest)**
```bash
# 1. Install VS Code
# 2. Install "Lean 4" extension
# 3. Extension will install Lean automatically
```

**Option 2: Command Line**
```bash
curl https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -sSf | sh
elan toolchain install stable
elan default stable
```

### Setup Your Project

```bash
cd /data1/yizhangh/ldp-pq
./setup_lean_project.sh
```

### Check Your File

**In VS Code:**
- Open `ButterflyCounting_2_2.lean`
- Lean checks automatically (green ✓ = good, red × = error)

**Command Line:**
```bash
lean --check ButterflyCounting_2_2.lean
```

## File Status Indicators

- ✅ **Green checkmark** = Proof complete and correct
- ❌ **Red X** = Error (syntax or type error)
- ⚠️ **Yellow warning** = Incomplete proof (has `sorry`)

## Common Commands

```bash
# Check a file
lean --check file.lean

# Build project
lake build

# Update dependencies
lake update

# Clean build
lake clean && lake build
```

## Your Current File

**File**: `ButterflyCounting_2_2.lean`
- ✅ Compiles (no syntax errors)
- ⚠️ Has `sorry` (incomplete proofs)
- 📝 Ready for proof development

**To complete**: Replace `sorry` with actual proofs

## What Happens When You Run Lean?

1. **Lean reads your `.lean` file**
2. **Checks syntax and types**
3. **Verifies proofs** (if complete)
4. **Reports errors or success**

**Lean is NOT like Python** - it doesn't "run" your code. It:
- ✅ Checks that your code is correct
- ✅ Verifies your proofs are valid
- ✅ Reports what's wrong (if anything)

## Example Workflow

```bash
# 1. Write/edit your Lean file
vim ButterflyCounting_2_2.lean

# 2. Check it
lean --check ButterflyCounting_2_2.lean

# 3. See errors (if any) and fix them

# 4. Repeat until it compiles

# 5. Complete proofs (replace sorry)
#    - Open in VS Code for best experience
#    - Lean will check in real-time
```

## Key Points

- **Lean files = Source code** (like `.py` or `.cpp` files)
- **You write them** (not generated)
- **Lean verifies them** (checks correctness)
- **VS Code is best** for interactive development
- **Command line works** for batch checking

## Need Help?

- **Documentation**: See `LEAN_SETUP_AND_USAGE.md`
- **Lean Manual**: https://leanprover.github.io/lean4/doc/
- **Community**: https://leanprover.zulipchat.com/
