# ORLQP Library

## Installation

### Build from source

```
git clone https://github.com/JTylerBoylan/orlqp
cd orlqp && mkdir build
cd build && cmake ..
cmake --build . --target install
```

### Add to CMake

```
find_package(orlqp REQUIRED)
...
target_link_libraries(my_executable PUBLIC orlqp::orlqp ...)
```