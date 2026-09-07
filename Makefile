FC = gfortran

FFLAGS_BASE = -ffree-form -fimplicit-none -Werror=implicit-interface
FFLAGS_DEBUG   = -g -O0 -Wall -Wextra -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
FFLAGS_RELEASE = -O3 -DNDEBUG
FFLAGS = $(FFLAGS_BASE) $(FFLAGS_RELEASE)

SRCDIR  = src
OBJDIR  = build/obj
APPDIR  = app
TESTDIR = test

OBJ_LIB = $(OBJDIR)/precision_mod.o \
          $(OBJDIR)/function_interfaces.o \
          $(OBJDIR)/test_functions.o \
          $(OBJDIR)/Gaussian_elimination.o \
          $(OBJDIR)/derivative_matrix.o \
          $(OBJDIR)/newton-nd.o

TARGET   = newton-nd
TEST_BIN = test_jacobian test_newtone

.PHONY: all clean main test debug release

all: main
main: $(TARGET)
test: $(TEST_BIN)

debug:
	$(MAKE) FFLAGS="$(FFLAGS_BASE) $(FFLAGS_DEBUG)"

release:
	$(MAKE) FFLAGS="$(FFLAGS_BASE) $(FFLAGS_RELEASE)"

$(TARGET): $(OBJ_LIB) $(OBJDIR)/main.o
	$(FC) $(FFLAGS) -o $@ $^

$(TEST_BIN): %: $(OBJ_LIB) $(OBJDIR)/%.o
	$(FC) $(FFLAGS) -o $@ $(OBJ_LIB) $(OBJDIR)/$*.o

$(OBJDIR)/precision_mod.o: $(SRCDIR)/precision_mod.f90 | $(OBJDIR)
$(OBJDIR)/function_interfaces.o: $(SRCDIR)/function_interfaces.f90 $(OBJDIR)/precision_mod.o | $(OBJDIR)
$(OBJDIR)/test_functions.o: $(SRCDIR)/test_functions.f90 $(OBJDIR)/precision_mod.o $(OBJDIR)/function_interfaces.o | $(OBJDIR)
$(OBJDIR)/Gaussian_elimination.o: $(SRCDIR)/Gaussian_elimination.f90 $(OBJDIR)/precision_mod.o | $(OBJDIR)
$(OBJDIR)/derivative_matrix.o: $(SRCDIR)/derivative_matrix.f90 $(OBJDIR)/precision_mod.o $(OBJDIR)/function_interfaces.o $(OBJDIR)/test_functions.o | $(OBJDIR)
$(OBJDIR)/newton-nd.o: $(SRCDIR)/newton-nd.f90 $(OBJDIR)/precision_mod.o $(OBJDIR)/derivative_matrix.o $(OBJDIR)/function_interfaces.o $(OBJDIR)/Gaussian_elimination.o | $(OBJDIR)

$(OBJDIR)/main.o: $(APPDIR)/main.f90 $(OBJDIR)/precision_mod.o $(OBJDIR)/newton-nd.o $(OBJDIR)/test_functions.o | $(OBJDIR)
$(OBJDIR)/test_jacobian.o: $(TESTDIR)/test_jacobian.f90 $(OBJDIR)/precision_mod.o $(OBJDIR)/derivative_matrix.o $(OBJDIR)/function_interfaces.o $(OBJDIR)/test_functions.o | $(OBJDIR)
$(OBJDIR)/test_newtone.o: $(TESTDIR)/test_newtone.f90 $(OBJDIR)/precision_mod.o $(OBJDIR)/newton-nd.o $(OBJDIR)/test_functions.o | $(OBJDIR)

$(OBJDIR)/%.o:
	$(FC) $(FFLAGS) -c -o $@ $<

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) $(TARGET) $(TEST_BIN) result.dat *.mod
