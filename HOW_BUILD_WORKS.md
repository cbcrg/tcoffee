# How T-Coffee build works

T-Coffee is automatically build by using the [Circle CI](https://circleci.com/gh/cbcrg/tcoffee) service.

The Circle build is controlled by the [circle.yml](.circleci/config.yml) file. 

## Build process main tasks 

1. The T-Coffee [build/build.sh](build/build.sh) script is executed in the [cbcrg/tcoffee-build-box:1.2](docker/Dockerfile.buildbox) container

2. The T-Coffee binaries produced by the build process are created in the folder 
  `~/publish/sandbox/build`. This folder is moved under path `build/tcoffee`
  
3. A new Docker image named `xcoffee` is created copying the content of the folder `build/tcoffee`. 
  The image is built by using this [docker/Dockerfile](docker/Dockerfile). 
  
4. The new `xcoffee` build is tested by running the [docker/run-tests.sh](docker/run-tests.sh) 
  tests suite. 
  
5. Test results are published to the following [location](http://www.tcoffee.org/Packages/tests/).

6. If all tests are green the `xcoffee` image is pushed to the [Docker Hub](https://hub.docker.com/r/cbcrg/tcoffee/tags/) 
  with the names `cbcrg/tcoffee:latest` and `cbcrg/tcoffee:<version.commit-id>`

7. The environment variable `RELEASE=0|1` is used to mark the build as beta or stable 
  (use the [build/make_release.sh] to trigger a new release build).
