PROG := sim.out
SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:%.cpp=%.o)
DEPS := $(SRCS:%.cpp=%.d)

CC := g++
CCFLAGS := -std=c++20
INCLUDEPATH := -I /usr/include/eigen3

all: $(DEPENDS) $(PROG)

$(PROG): $(OBJS)
	$(CC) $(CCFLAGS) -o $@ $^ -O3 -fopenmp -mavx

.cpp.o:
	$(CC) $(CCFLAGS) $(INCLUDEPATH) -MMD -MP -MF $(<:%.cpp=%.d) -c $< -o $(<:%.cpp=%.o) -O3 -fopenmp -mavx

.PHONY: clean
clean:
	$(RM) $(PROG) $(OBJS) $(DEPS)

-include $(DEPS)