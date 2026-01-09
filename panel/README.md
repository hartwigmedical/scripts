# Publish new preprocess-recon-cnv release

First make sure you can publish to the Hartwig Docker repository:

    gcloud auth configure-docker europe-west4-docker.pkg.dev

Build the local Docker image:

    docker build -t preprocess-recon-cnv .

Test the local Docker image:

    docker run preprocess-recon-cnv

If the local Docker image works as expected, tag it with the public name 
and a [semantic version](https://semver.org) (we use `0.0.1` as the example version):

    docker tag preprocess-recon-cnv europe-west4-docker.pkg.dev/hmf-build/hmf-docker/preprocess-recon-cnv:0.0.1

Finally, push the Docker image to the Hartwig Docker repository:

    docker push europe-west4-docker.pkg.dev/hmf-build/hmf-docker/preprocess-recon-cnv:0.0.1
