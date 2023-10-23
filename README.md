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

## Installation

### Build from source

```
git clone https://github.com/JTylerBoylan/orlqp
cd orlqp && mkdir build
cd build && cmake ..
make install
```

### Add to CMake

```
find_package(orlqp REQUIRED)
...
target_link_libraries(my_executable PUBLIC orlqp::orlqp ...)
```

## Docker (Dev)

```
sudo chmod +x run_dev.sh
./run_dev.sh
```

*Or `./cuda_run_dev.sh` for large matrices (thousands of decision variables)*