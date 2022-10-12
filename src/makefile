PROG := sim.out
SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:%.cpp=%.o)
DEPS := $(SRCS:%.cpp=%.d)
OMPS := -O3 -fopenmp -mavx

CC := g++
CCFLAGS := -std=c++17
INCLUDEPATH := -I /usr/include/eigen3 -I /usr/local/include/eigen3
DBG := -g3 -Wall

all: $(DEPENDS) $(PROG)

$(PROG): $(OBJS)
	$(CC) $(DBG) $(CCFLAGS) -o $@ $^ $(OMPS)

.cpp.o:
	$(CC) $(DBG) $(CCFLAGS) $(INCLUDEPATH) -MMD -MP -MF $(<:%.cpp=%.d) -c $< -o $(<:%.cpp=%.o) $(OMPS)

.PHONY: clean
clean:
	$(RM) $(PROG) $(OBJS) $(DEPS)

-include $(DEPS)
