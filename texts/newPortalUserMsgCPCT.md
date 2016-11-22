# Hartwig Medical Download Portal

Dear client,

We are about to setup an account for you at our download portal. This email provides you with information to have a quick start.

Regards,  
The HMF Team

### DISCLAIMER CPCT-02

Data collected as part of the CPCT-02 protocol is only to be used to investigate specific research questions within the scope of this protocol and according to the study contract.

Distribution and/or sharing of these data is not allowed without prior approval of the sponsor. Approval from the sponsor is required prior to distribution of data to countries outside of the European Union (also for anonymized or coded data).

Please be aware that these data are for internal use only.

**Scope of the CPCT-02 protocol**
The primary objective of this study is to analyze the individual cancer genome in cancer patients to develop future predictors for response to systemic treatment in individual patients with cancer.

Secondary objectives of this study are:
- To determine the amount of biopsy samples with sufficient DNA for analysis
- To determine the amount of biopsy samples with an adequate mutational profile
- To collect and anonymously interpret all mutational profiles obtained using this protocol
- To determine changes in the mutational profile under the influence of systemic treatment
- To explore and analyze the individual microRNA, (phospho)proteomic profiles and organoid cultures in patients with cancer to develop future predictors for response to systemic treatment
- To explore the correlation between mutational profiles in solid tumor biopsies and liquid biopsies (circulating tumor DNA)


### Steps to get access to the HMF portal
1. You will receive two emails from Schuberg Philis containing:
    - your username
    - link to start the token enrollment process
2. Follow the instructions in the enrollment email carefully. You will need to open the link in the email on your smartphone to install an app called "MobilePass". In the app you will choose a PIN-code that you will use in the future to generate a passcode in this app. Remember this PIN-code. With your username and a generated passcode you can login to the portal.
3. When the setup is done, login at: https://connect.schubergphilis.com/hmfdownloads


### What to do after successful login
Once logged in you will have access to data from your group or center. For each sample/patient certain files can be downloaded by creating links (valid for 24 hours). The links can be used to download (via browser or standard tools like WGET/cURL) or directly read into a program like IGV (Integrative Genome Browser).

**Example download using WGET:**
- create the links (by clicking the most right icon next to a run)
- copy the links into a new text file (eg download.txt)
- use the following command to download: 
```sh
wget --content-disposition -i download.txt
```

**Example loading BAM file in IGV:**
- create the links for the BAM and accompanying BAI file (by clicking the most right icon next to a run)
- open the IGV program and choose "Load from URL"
- paste the BAM file link at "File URL"
- paste the BAI file link "Index URL"

For more information and contact visit our website:
[Hartwig Medical Foundation](https://www.hartwigmedicalfoundation.nl)
