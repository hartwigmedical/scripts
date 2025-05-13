## Build Docker image
First, check Artifact Registry in hmf-build for what version numbers have already been used for that Docker image. Pick an unused version number.

Build an image and push it:
```shell script
docker build . -t europe-west4-docker.pkg.dev/hmf-build/hmf-docker-crunch/${name}:${version}
docker push europe-west4-docker.pkg.dev/hmf-build/hmf-docker-crunch/${name}:${version}
```
e.g:
```shell script
docker build . -t europe-west4-docker.pkg.dev/hmf-build/hmf-docker-crunch/auto-compare:0.0.9
docker push europe-west4-docker.pkg.dev/hmf-build/hmf-docker-crunch/auto-compare:0.0.9
```
Use hmf-docker instead of hmf-docker-crunch for more official workflows, but only when first discussed with the INFRA team.

## Push new workflow
Check whether there is already a version of this workflow pushed through http://hawe.verification-1/swagger-ui/index.html or
```shell script
curl -X 'GET' \
  'http://hawe.verification-1/workflows' \
  -H 'accept: */*'
```

Make sure to uptick the version before pushing.

Push:
```shell script
curl -X PUT -H "Content-Type: application/json" \
  --data-binary "@${workflow}" \
  'http://hawe.verification-1/workflow'
```
If this fails, possibly first delete current version of the workflow from hawe after discussing this with the people that might be using it.