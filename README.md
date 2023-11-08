# ORLQP Library

## Dependencies

### CMake

[Website](https://cmake.org/)

```
sudo apt-get install -y cmake
```

### Eigen

[Website](https://eigen.tuxfamily.org/index.php?title=Main_Page)

```
sudo apt-get install -y libeigen3-dev
```

### OSQP

[Github](https://github.com/osqp/osqp)

[Website](https://osqp.org/)

### GiNaC

[Website](https://www.ginac.de/)

```
sudo apt-get install -y libginac-dev libcln-dev
```

## Installation

### Build from source

```
git clone https://github.com/JTylerBoylan/orlqp
cd orlqp && mkdir build
cd build && cmake ..
make install
```

### Build in Docker

```
# Check for updates
ADD https://api.github.com/repos/JTylerBoylan/orlqp/git/refs/heads/main version.json

# Clone and compile library
RUN git clone https://github.com/JTylerBoylan/orlqp
RUN cd orlqp && mkdir build
RUN cd orlqp/build && cmake ..
RUN cd orlqp/build && cmake --build . --target install
```

*NOTE: You may need to add `RUN ldconfig` to the Dockerfile if you get a linking error.*

## CMake

```
find_package(orlqp REQUIRED)
...
target_link_libraries(my_executable PUBLIC orlqp::orlqp ...)
```