# PR previews

For setting up the pr previews we need to have AWS and the github repository configured.

## AWS 
Create a new IAM user with the following permissions:

```
  {
	"Version": "2012-10-17",
	"Statement": [
            {
                "Sid": "Statement1",
                "Effect": "Allow",
                "Action": [
                    "s3:CreateBucket",
                    "s3:DeleteBucket",
                    "s3:GetObject",
                    "s3:PutObject",
                    "s3:DeleteObject",
                    "s3:PutBucketWebsite",
                    "s3:PutObjectAcl",
                    "s3:ListBucket",
                    "s3:PutBucketPublicAccessBlock",
                    "s3:PutBucketOwnershipControls",
                    "s3:ListMultipartUploadParts",
                    "s3:AbortMultipartUpload",
                    "s3:GetObjectAcl"
                ],
                "Resource": [
                    "arn:aws:s3:::pr-preview-bioconductor-*"
                ]
            }
	    ]
    }
```

Make sure that you keep a hold of the `Access Key ID` and the `Secret Access Key` for that IAM user, we  will need them to set up github.


## Github Actions

There are two yaml files in the workflows directory associated with the pr previews which are:

`pr_deploy.yaml` and `pr_close.yaml`

The former will run the `Pr - Preview` action on any non draft pull request, and the latter will run the job `PR - Closed` action on any closed or non-open pull request. Below are the jobs and steps for each action 

### pr_deploy.yaml *(PR - Preview)*
#### dev-pr-create-s3:
##### steps: 
- Checks out repository.
- Caches and retrieves the output directory(for faster build times).
- Sets up and builds the website to the output directory.
- Creates bucket, if bucket does not already exist.
    - This step was required as the next action had outdated methods for creating the bucket.
- Deploys contents of the output directory to bucket.

##### Repository configuration for this job:

- Go to the settings page for the repository
    - In the repository page click on the 'Settings' tab at the top of the page, under the repository name.
- Create a new environment called "dev".
    - On the left side panel, click on environment.
    - There should be a button saying "new environment" on the page, once clicked we enter in our environment name `dev` and click "configure environment". *( This is to make sure that the access keys, which will be added in the next step are only accessible on a dev environment.)*
- Add your IAM access key and ID *(from the aws section)* as secrets to the dev environment.
    - Once on the environments page there should be a section called "environment secrets". 
    - Click on "add secrets", give the secret a name of "AWS_SECRET_ACCESS_KEY" and paste in the IAM access key into the 'value' field *(from the aws section)*.
    - Once again click "add secrets" and give this secret a name of "AWS_ACCESS_KEY_ID" and paste in the IAM access key ID into the 'value' field *(from the aws section)*.

#### kpi_metrics:

##### steps: 
- Checks out repository.
- Installs node for running the KPI scripts, and the relevant packages (`psi` and `@axe-core/cli`)
- Runs PageSpeed insights:
    - Using node, we run the js file PageSpeed.js and provide it with the relevant parameters:
        - `url` - website url 
        - `strategy` - desktop or mobile 
        - `threshold` - range from 0 - 100, the minimum score you would like to pass.
        - `apiKey` - api key for pageSpeed this can be "none" or your api key if you have one.
    - This step is run for both mobile and desktop, the output is stored as an environment variable as `DESKTOP_PSI_RESULTS` and `MOBILE_PSI_RESULTS` respectively.
- Add pageSpeed comments to the pull request.
    - This step gets the outputs of the pageSpeed step from the github environment variable and comments them onto the current pull request. We provide a `message`, this is the content we what to comment, and a `message-id`, since we will also make a comment for the AXE results.
- Run Axe accessibility:
    - Using npx, we run the command 'npx run axe (url)' using the bioconductor url, the output is stored in the github environment variable `AXE_RESULTS`.
- Add AXE comments to pull request.
    - Similar to the pageSpeed comments step, we comment Axe results onto the pull request, the `message` this time being the contents of the `AXE_RESULTS` env variable. We provide also different `message-id`.

*There is no specific repository configuration for this job.*

### pr_close.yaml *(PR - Closed)*
##### steps:
- You need the secrets thats been used to configure the dev environment for pr_deploy
- Checks the bucket to see if it exists, to do that you need:
    - access key ID for iam user
    - access key for iam user
- If it does exist, an environment variable BUCKET_EXISTS = true is stored 
- Check to see if the bucket exists in AWS, to do that you need:
    - github token for the current PR
- If the bucket does exist, then it will be deleted from AWS


*For repository configuration, this job needs the environment `dev` and the AWS secrets we setup for that environment earlier, therefore there will no other setup required.*
