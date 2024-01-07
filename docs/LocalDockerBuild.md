### Local ldmxsw docker build and push

You can make changes to and build the base ldmxsw container (the one used for developing and running ldmx-sw locally). To do so:
1. Pull down the ldmxsw docker repo: `git clone https://github.com/LDMX-Software/docker.git`
2. Move into that directory (`cd docker`, or whatever you called it).
3. Make the desired changed, typically to the `Dockerfile` file.
4. Build using docker with the command `docker build --build-arg NPROC=<nproc> -t <tag> ./` where <tag> is the name of the tag that you would like to give it and `<nproc>` is the number of parallelized processors you want to use in the build to make it go faster (this is optional: defaults to 1).

For the tag, to use this container locally, you will likely want to name it `lmdx/local:<tag>`, such that when setting up the ldmxsw environment, you can do `ldmx use local <tag>`.

You can always create a new tag name by doing `docker tag SOURCE_IMAGE[:TAG] TARGET_IMAGE[:TAG]`.

You can push your build to DockerHub or to other image repository doing the usual `docker push <tag>`. E.g. for pushing to Github IO, you can do something like:
```
#make the tag
docker tag ldmx/local:latest ghcr.io/<gh_username>/ldmxsw:latest

#push
docker push ghcr.io/<gh_username>/ldmxsw:latest
```
