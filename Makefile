CXX = clang++
CXXFLAGS += -std=c++17 -O0 -g
CXXFLAGS += -Werror -Wall -Wextra -Wshadow -Wpedantic
CXXFLAGS += -Wno-unused-function -Wno-unused-parameter -Wno-unused-variable # TODO: disable for production builds

OBJECTS := main.o common.o matrix.o transform.o camera.o film.o geometry.o light.o spectrum.o sampling.o bsdf.o intersection.o entity.o scene.o
OBJECTS := $(addprefix build/, $(OBJECTS))

.PHONY: default
default: dirs pt

dirs:
	@mkdir -p bin build .deps

pt: $(OBJECTS)
	@$(CXX) $(LDFLAGS) $^ -o bin/pt

.SECONDARY:
build/%.o: src/%.cpp
	@$(CXX) $(CXXFLAGS) -o $@ -c $< -MMD -MP -MF .deps/$*.d

clean:
	@rm -rf bin build .deps

-include .deps/*.d
