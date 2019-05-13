#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------
EXENAME := driver.exe

NETCDF := /usr/local/netcdf-4.7.0

RM := rm
FC := pgf90

SRCLIST := Base_Srclist_f
dirs := .:../share:../rad_xfer:../photo_rates:../molec_ox_xsect:../../MusicBox/MusicBox_host/src

# Expand any tildes in directory names. Change spaces to colons.
VPATH    := $(foreach dir,$(dirs),$(dir)) 
VPATH    := $(subst $(space),:,$(VPATH))               

# Get list of files and determine objects and dependency files
base_srclist_f := $(shell cat $(SRCLIST))
OBJS      := $(foreach file,$(base_srclist_f),$(file:.F90=.o))

FFLAGS := -c -g -Ktrap=fp -Mbounds -Kieee -traceback -I$(NETCDF)/include 
LDFLAGS := -L$(NETCDF)/lib -lnetcdff -L$(NETCDF)/lib -lnetcdf -Wl,-rpath -Wl,$(NETCDF)/lib

all: $(EXENAME)

#------------------------------------------------------------------------
# Default rules
#------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .F90 .o

.F90.o:
	$(FC) $(FFLAGS) -o $@ $<

$(EXENAME):	$(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) -f *.o *.mod $(EXENAME) *.out TUV.diags fort.10

