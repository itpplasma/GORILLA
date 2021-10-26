LATEST_PFUNIT_DIR := $(PFUNIT_DIR)/
include $(LATEST_PFUNIT_DIR)/include/PFUNIT.mk

all: test_various_functions

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I../../OBJS

test_various_functions_TESTS := test_various_functions_mod.pf
test_various_functions_OTHER_LIBRARIES :=  -L../../OBJS -lsut
$(eval $(call make_pfunit_test,test_various_functions))


clean:
	$(RM) *.o *.mod *.a  *.inc
	$(RM) test_various_functions_mod.F90 test_various_functions

