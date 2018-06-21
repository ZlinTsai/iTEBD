CONFIG_FILE := make.inc
SRCINFO     := srcinfo.inc
ifeq ($(wildcard $(CONFIG_FILE)),)
    $(error $(CONFIG_FILE) not found.)
endif

include $(CONFIG_FILE)
include $(SRCINFO)

# The name of binary object.
EXU := $(MAINSCRIPT)

# The file names of source codes.
src := $(SRC)

COMMON_FLAGS_INCS += $(foreach root,$(ROOTS),-I$(root)/include)
COMMON_FLAGS_LIBS += $(foreach root,$(ROOTS),-L$(root)/lib)
CXXFLAGS += $(UNI10CXXFLAGS)
UNI10LDLIBRARY := -luni10

all: job.exe 

job.exe: $(EXU)
	$(CXX) $(COMMON_FLAGS_INCS) $(COMMON_FLAGS_LIBS) $(CXXFLAGS) $^ -o $@ $(UNI10LDLIBRARY) 

clean:
	rm -rf test.exe
