######################################################################
LIBMESH_DIR ?= /Users/vasvav/Work/libmesh

include $(LIBMESH_DIR)/Make.common

target   := ./rdcFEs.$(METHOD)
srcfiles := $(wildcard src/*.C)
objects  := $(patsubst src/%.C, src/%.$(obj-suffix), $(srcfiles))
###############################################################################

.PHONY: dust clean distclean

###############################################################################

all:: $(notdir $(target))

# Production rules:  how to make the target - depends on library configuration
$(notdir $(target)): $(objects)
	@echo "Linking "$@"..."
	@$(libmesh_LIBTOOL) --tag=CXX $(LIBTOOLFLAGS) --mode=link \
	  $(libmesh_CXX) $(libmesh_CXXFLAGS) $(objects) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS) $(EXTERNAL_FLAGS)

###############################################################################
dust:
	@echo "Deleting old output and runtime files"
	@rm -f out*.m job_output.txt output.txt* *.gmv.* *.plt.* out*.xdr* out*.xda* PI* *.e *.ex2 *.vtk *.vtu *.pvd *.pvsm $(target)

clean: dust
	@rm -f $(objects) *.$(obj-suffix) *.lo

clobber: clean 
	@rm -f $(target)

distclean: clean
	@rm -rf *.o .libs

echo:
	@echo srcfiles = $(srcfiles)
	@echo objects = $(objects)
	@echo target = $(target)

run: complete

complete: $(wildcard *.in)
#	@$(MAKE) dust
	@$(MAKE) -C $(dir $(target)) $(notdir $(target))
	@echo "***************************************************************"
	@echo "* Running App " $(notdir $(target))
	@echo "***************************************************************"
	@echo " "
	mpiexec -n 1 $(target) ${LIBMESH_OPTIONS} 2>&1 | tee output.txt
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running App " $(notdir $(target))
	@echo "***************************************************************"

###############################################################################
