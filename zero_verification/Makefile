CXX      := g++
CXXFLAGS := -Wall -Wextra -std=c++17

TARGETS  := 41_upper 41_lower

all: $(TARGETS)

41_upper: 41_upper.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

41_lower: 41_lower.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(TARGETS)

.PHONY: all clean
