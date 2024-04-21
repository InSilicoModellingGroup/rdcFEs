######################################################################
LIBMESH_DIR = /Users/vasvav/Public/libmesh-d3bda6c7009c6b3241ef2e6c999eff577116dc68

include $(LIBMESH_DIR)/Make.common

this_project = $(shell pwd)
target   := $(this_project)/./rdcFEs.$(METHOD)
srcfiles := $(wildcard $(this_project)/src/*.C)
objects  := $(patsubst %.C, %.$(obj-suffix), $(srcfiles))
###############################################################################

.PHONY: dust clean distclean

###############################################################################

all: $(notdir $(target))

# Production rules:  how to make the target - depends on library configuration
$(notdir $(target)): $(objects)
	@echo "Linking "$@"..."
	@$(libmesh_LIBTOOL) --tag=CXX $(LIBTOOLFLAGS) --mode=link \
	  $(libmesh_CXX) $(libmesh_CXXFLAGS) $(objects) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS) $(EXTERNAL_FLAGS)

###############################################################################
dust:
	@echo "Deleting old output and runtime files"
	@rm $(this_project)/src/*.o
	@rm -f out*.m job_output.txt output.txt* *.gmv.* *.plt.* out*.xdr* out*.xda* PI* *.e *.ex2 *.vtk *.vtu *.pvd $(target)

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

###############################################################################
