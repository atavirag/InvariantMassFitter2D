# Makefile for InvariantMassFitter with organized build structure

CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2 -fPIC
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs) -lRooFit -lRooFitCore -lMinuit

# Directory structure
SRCDIR = src
BUILDDIR = build
DATADIR = data
PLOTDIR = plots

# Include directories
INCLUDES = -I$(SRCDIR) $(ROOTFLAGS)

# Source files
SOURCES = $(SRCDIR)/InvariantMassFitter.cxx
HEADERS = $(SRCDIR)/InvariantMassFitter.h
OBJECTS = $(BUILDDIR)/InvariantMassFitter.o

# Executables
TARGET = $(BUILDDIR)/runFit

# Create directories if they don't exist
$(shell mkdir -p $(BUILDDIR) $(DATADIR) $(PLOTDIR))

# Default target
all: $(TARGET)

# Compile object files
$(BUILDDIR)/%.o: $(SRCDIR)/%.cxx $(HEADERS)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Link executable
$(TARGET): $(OBJECTS) $(SRCDIR)/runFit.cxx
	@echo "Linking $@..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJECTS) $(SRCDIR)/runFit.cxx -o $@ $(ROOTLIBS)
	@echo "Build complete! Executable: $(TARGET)"

# Compile as shared library (optional)
lib: $(OBJECTS)
	@echo "Creating shared library..."
	$(CXX) -shared -o $(BUILDDIR)/libInvariantMassFitter.so $(OBJECTS) $(ROOTLIBS)

# Run the analysis
run: $(TARGET)
	@echo "Running analysis..."
	cd $(DATADIR) && ../$(TARGET)

# Clean build files only
clean:
	@echo "Cleaning build artifacts..."
	rm -f $(BUILDDIR)/*.o $(BUILDDIR)/*.so $(BUILDDIR)/*.d $(BUILDDIR)/runFit
	@echo "Build directory cleaned"

# Clean everything including output
cleanall: clean
	@echo "Cleaning all generated files..."
	rm -f $(DATADIR)/*.root
	rm -f $(PLOTDIR)/*.pdf $(PLOTDIR)/*.png
	@echo "All generated files removed"

# Clean only data files
cleandata:
	@echo "Cleaning data files..."
	rm -f $(DATADIR)/*.root

# Clean only plots
cleanplots:
	@echo "Cleaning plots..."
	rm -f $(PLOTDIR)/*.pdf $(PLOTDIR)/*.png

# Show build info
info:
	@echo "Build Configuration:"
	@echo "  Source dir:  $(SRCDIR)"
	@echo "  Build dir:   $(BUILDDIR)"
	@echo "  Data dir:    $(DATADIR)"
	@echo "  Plot dir:    $(PLOTDIR)"
	@echo "  CXX:         $(CXX)"
	@echo "  ROOT:        $(shell root-config --version)"

# Help
help:
	@echo "Available targets:"
	@echo "  all       - Build executable (default)"
	@echo "  lib       - Build shared library"
	@echo "  run       - Build and run analysis"
	@echo "  clean     - Remove build files"
	@echo "  cleanall  - Remove build files and output"
	@echo "  cleandata - Remove only data files"
	@echo "  cleanplots- Remove only plot files"
	@echo "  info      - Show build configuration"
	@echo "  help      - Show this message"

.PHONY: all clean cleanall cleandata cleanplots lib run info help