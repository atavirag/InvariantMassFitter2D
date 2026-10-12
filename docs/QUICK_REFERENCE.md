# Quick Reference Card

## 📁 File Organization Cheat Sheet

```
✅ COMMIT TO GIT          ❌ DON'T COMMIT (build/)      ❌ DON'T COMMIT (data/)
src/*.{h,cxx}             build/*.so                    data/*.root
python/*.py               build/*.d
Makefile                  build/*.o                     ❌ DON'T COMMIT (plots/)
.gitignore                build/runFit                  plots/*.pdf
README.md                 *_cxx.so                      plots/*.png
                          *_cxx.d
```

## ⚡ Common Commands

### Building
```bash
make            # Compile
make clean      # Clean build files
make run        # Build and execute
make cleanall   # Nuclear option
```

### Running
```bash
./build/runFit                                    # C++ fitting
python python/visualize_with_fits.py              # Python plots
python python/check_analysis.py --all            # Diagnostics
scripts/run_complete_analysis.sh                  # Full pipeline
```

### Cleaning
```bash
make clean                # Build files only
make cleandata           # ROOT files only
make cleanplots          # Plots only
make cleanall            # Everything generated
```

## 🔍 Where Are My Files?

| Looking for... | Check... |
|----------------|----------|
| Source code | `src/` |
| Compiled libraries | `build/*.so` |
| Executables | `build/runFit` |
| ROOT files | `data/*.root` |
| Plots | `plots/*.pdf` |
| Python scripts | `python/*.py` |
| Documentation | `docs/` |

## 🎯 File Extensions Guide

| Extension | What is it? | Keep in Git? | Where? |
|-----------|-------------|--------------|--------|
| `.h` | C++ header | ✅ Yes | `src/` |
| `.cxx` | C++ source | ✅ Yes | `src/` |
| `.o` | Object file | ❌ No | `build/` |
| `.so` | Shared library | ❌ No | `build/` |
| `.d` | Dependencies | ❌ No | `build/` |
| `.root` | ROOT data | ❌ No | `data/` |
| `.pdf` | Plot | ❌ No | `plots/` |
| `.py` | Python script | ✅ Yes | `python/` |
| `.sh` | Shell script | ✅ Yes | `scripts/` |
| `.md` | Documentation | ✅ Yes | `docs/` |

## 🚨 Common Mistakes

### ❌ DON'T DO THIS
```bash
git add *.so          # Never commit .so files!
git add data/         # Never commit large data!
git add build/        # Never commit build artifacts!
```

### ✅ DO THIS INSTEAD
```bash
git add src/          # Source code
git add python/       # Scripts
git add Makefile      # Build config
git status --ignored  # Check what's ignored
```

## 🛠️ First Time Setup

```bash
# 1. Organize files
./setup_project.sh

# 2. Initialize git
git init
git add .
git commit -m "Initial commit"

# 3. Build
make

# 4. Run
make run
```

## 🔄 Daily Workflow

```bash
# Morning
git pull
make clean
make

# Work
vim src/MyFile.cxx
make
./build/runFit

# Check
python python/check_analysis.py

# Commit
git add src/MyFile.cxx
git commit -m "Added feature X"
git push
```

## 📦 When Sharing Code

### What to Share
```bash
# Create clean archive
git archive --format=tar.gz --output=analysis.tar.gz HEAD

# Or zip source only
zip -r analysis_src.zip src/ python/ Makefile README.md .gitignore
```

### What NOT to Share
- `build/` directory
- `data/*.root` files (send separately if needed)
- `plots/` (they'll regenerate)
- Any `*.so`, `*.o`, `*.d` files

## 🐛 Quick Fixes

### "Cannot find file"
```bash
# Check you're in project root
pwd

# Check directory structure
ls -R | grep -E "^(src|build|data):"
```

### "Build fails"
```bash
make clean
make info   # Check configuration
make        # Rebuild
```

### "Git tracking wrong files"
```bash
# Remove from git (keep on disk)
git rm --cached build/*.so
git rm --cached data/*.root

# Update .gitignore
echo "build/" >> .gitignore
echo "data/" >> .gitignore

git add .gitignore
git commit -m "Fix gitignore"
```

## 📝 Quick Edits

### Change where files go
```makefile
# In Makefile, change:
BUILDDIR = build          # → your_build_folder
DATADIR = data            # → your_data_folder
PLOTDIR = plots           # → your_plots_folder
```

### Add file to git after ignoring
```bash
# In .gitignore, add exception:
!data/important_test.root

# Then:
git add data/important_test.root
```

## 💾 Disk Space Management

```bash
# Check sizes
du -sh build/ data/ plots/

# Free up space
make cleanall              # Remove all generated files
rm -rf build/ data/ plots/ # Nuclear option

# Rebuild when needed
make
```

## 🎓 Learning More

| Topic | File |
|-------|------|
| C++ Fitting | `docs/README_Fitter.md` |
| Python Viz | `docs/README_Python_Viz.md` |
| Full Guide | `docs/PROJECT_OVERVIEW.md` |
| This File | `QUICK_REFERENCE.md` |

## 🎯 Goals

- ✅ Clean separation of source vs generated
- ✅ Easy to build: `make`
- ✅ Easy to clean: `make clean`
- ✅ Git-friendly: Only source in repo
- ✅ Reproducible: Anyone can rebuild

---

**Print this and keep it handy!** 📋