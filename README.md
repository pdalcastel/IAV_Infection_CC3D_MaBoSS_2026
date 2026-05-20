# Development of a digital innate immune twin for influenza A infection of the airway

This repository contains the source code, results, raw data and other scripts of the IAV-CC3D-MaBoSS model and the fitting of this model to experimental data. This model is an agent-based virtual tissue representation of the lung epithelium with virus-host antagonism pathways, semi-infectious viral particles and genome complementation on superinfection. The manuscript describes the model, full citation below,

> Waiting for review

---

## Model Description

There are key biological observations incompatible with current mechanistic models of IAV infection, for example, the robust IFN response of epithelial cells at low MOI regimes despite the observed effectiveness of NS1 - determined in high MOI experiments. To accurately reproduce IAV infection behavior and its effect on IFN response for simulating prevention and treatment, we decided to incorporate biological hypothesis from the literature that no other mechanistic models, known to the authors, have considered. The main hypothesis that drastically distinguish this model from its predecessors are 1) viruses carry gene defects, 2) once inside the cell, different viral genomes complement each other, 3) viruses that are NS1 competent inhibit IFN response. With these hypothesis, our model could reproduce both the observed effectiveness of NS1 at high MOI, and the robust IFN response at low MOI. These hypothesis require that the model is agent-based, because each cell tracks the gene defects of each virus it endocytoses. Space is necessary to represent the spread of infection and IFN response, which are both localized and depend on whether cells are infected by NS1 incompetent viruses or not. For more details, check out the manuscript.

The model structure considers:
* **Universe**: 2D square lattice
* **Cell types**: uninfected, infected, dead
* **Diffusive fields**: Virus field, IFN field
* **Signaling**: intracellular network model, conditional secretion of IFN, conditional endocytose of viruses

See full model structure in the manuscript.

---

## Model behaviors

* **High MOI suppression of IFN response**: at high MOI, most cells are infected by NS1 competent virions, and therefore cannot robustly undergo IFN response
* **Low MOI robust IFN response**: because viruses can be defective, in a low MOI setting, cells are more likely to be infected by NS1 incompetent virions, becoming important sources of IFN production. Via the IRF7 pathway, IFN signaling is amplified and the cells undergo robust IFN response 
* **Virus endocytosis aggravates infection**: increasing virus endocytosis increases infection severity. Predecessor models predict the other way around. This is results from the virus reliance on complementatino for productive infection
* **IFN treatment**: increasing the concentration of IFN field prior to infection by providing IFN to the culture and increasing baseline IFN secretion increases antiviral behavior of cells, slowing down infection spread and decreasing virus production 

---

## Data and Code Availability

All model code, simulation files, analysis scripts, and data are openly available in this GitHub repository.

* **GitHub Repository (Development Version)**: up to date repository with latest changes
* **Zenodo Archive (Permanent Version)**: corresponds to the GitHub repository at the moment of paper submission 

---

## Running the Simulation Locally (installation of Conda and CompuCell3D required)

### 1. Setup with CompuCell3D Binaries
1. Follow installation instruction from the [Official Website](https://compucell3d.org/SrcBin) to download CompuCell3D 4.7.0
2. Once CompuCell3D is installed, launch CompuCell3D Player
3. Download the simulation folder `VRPDM_Continuous_Activation_ISG` from this repository by cloning the repository or by directly downloading the folder alone using the [https://download-directory.github.io/] app.
4. In the player, go to `File -> open simulation`, then navigate to the simulation folder `VRPDM_Continuous_Activation_ISG`, and open the `VRPDM_Continuous_Activation_ISG.cc3d` file.
6. Click the Play button

### 2. Environment Setup with Conda

Alternatively, you can obtain CompuCell3D via Conda. First, download [Anaconda](https://www.anaconda.com/products/distribution) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

Once Conda is installed, open the Anaconda Prompt (on Windows) or your terminal (on macOS/Linux) and run the following commands step-by-step:

```bash
# 1. Create a new conda environment named 'IAV_CC3D_MaBoSS' with Python 3.12. You can choose other environment name.
conda create -n IAV_CC3D_MaBoSS python=3.12

# 2. Activate the environment
conda activate IAV_CC3D_MaBoSS

# 3. Install Mamba (recommended)
conda install -c conda-forge mamba

# 4. Install CompuCell3D
mamba install -c main -c conda-forge -c compucell3d compucell3d=4.7.0

# 5. Open CC3D Player
python cc3d.player5

# 6. In the player, go to File -> open simulation, then navigate to the simulation folder "VRPDM_Continuous_Activation_ISG", and open the "VRPDM_Continuous_Activation_ISG.cc3d" file

# 7. Click the play button
```

## (NOT AVAILABLE) - Running the Simulation Online (no installations required)

The easiest way to run the simulation is by launching the simulation online in myBinder.org. 

1. Go to the folder `Jupyter_Notebook_Simulation` in this repository
2. Copy the browser link, or just copy this link: [https://github.com/pdalcastel/IAV_Infection_CC3D_MaBoSS_2025/tree/main/Jupyter_Notebook_Simulation]
3. Open [Binder](myBinder.org) in your navigator, myBinder.org
4. Paste the link on the field `GitHub repository name or URL` in the Binder website
5. Click `launch`
6. Once Jupyter opens, you can open the `JN_v4.ipynb` file from the project pannel
7. Run the `JN_v4.ipynb` file by clikcing the `Restart Kernel and Run All` button, which has a double right-arrow symbol. Alternatively, you can run cell by cell

If you running the JN_v4.ipynb cell by cell, you can modify the simulation parameters in the table, then click the button "Convert to Dictionary" to apply the changes.

---

## Replicating the Manuscript Results

#### Fitting procedure: 

1.  For the fitting, we set initial values of the parameters `secretion_rate_V=71.6/6.0`, `secretion_rate_IFN=2.5`, and `MOI_PFU=0.26`
2.  We run 5 replicates using CompuCell3D Player. The 5 replicates go in the output directory
3.  To get the averages and error bars, we use the avg.py script, which you can run using Python. You need to specify the name of the file you want to average without including the replicate number
4.  Open the OriginLab project in `Results_and_Plots -> final_calibration_with_correction.opju` (OriginLab required)
5.  Transfer the averaged file to the origin lab project. Plot virus vs time and IFN vs time together with the experimental data.
6.  Compare the simulation result to the experimental data. If the time series are to the left of the experimental curves, lower MOI. If the time series are to the right of the experimental curves, increase MOI. To raise/decrease IFN or virus levels, change the secretion of IFN and virus accordingly
7.  Iterate through steps 2-6, or simply use the fit parameters `secretion_rate_V=1900./6.`, `secretion_rate_IFN=4.55`, and `MOI_PFU=0.0012`
Of course this is an initial guess for the fitting parameter set, not necessarily the optimal fit. To assess the quality of the fit, we need to run a likelihood profile, a.k.a. a confidence region plot around the fit parameter set:
9.  Access a Slurm HPC cluster
10.  Transfer the cluster scripts in `cluster -> likelihood_profile`  to your cluster folder in the HPC
11.  You may need to set up Anaconda/Miniconda in your HPC cluster
12.  Install CompuCell3D in your cluster using the conda installation (see previous section or follow the steps in the [Official Website](https://compucell3d.org/SrcBin))
14.  Activate your CC3D environment
15.  Adapt the scripts to your system, such as paths, user name etc:
```bash
# Open the likelihood profile script

# Set the wanted number of replicas, total number of parameter sets, or leave as is
replicas = 10
N = 1000

# Set the maximum number of jobs. This depends on the specifics of your cluster
MAX_JOBS = 9000

# Set up and down limits for the three fitting parameters, or leave as is
parameters = {
'MOI_PFU': (p.MOI_PFU/5., p.MOI_PFU*5), 
'secretion_rate_V': (p.secretion_rate_V*0.5, p.secretion_rate_V*2.),
'secretion_rate_IFN': (p.secretion_rate_IFN/1.5, p.secretion_rate_IFN*1.5)
}
```
17. In the generate_sbatch_string function, change output, email and path to script
18. Run the script by calling `python likelihood_profile_script.py`
19. Collect the output files and run the `avg_likelihood_profile.py` script on the output. This script is inside the `Outpu.zip` file
20. Run the `error_calculation.py` script on the averaged files. You need to also have the `experimental_data.csv` file present
21. Open the `likelihood_profile.opju` project (OriginLab required) from the `Results_and_Plots` folder
22. Transfer the error table to the `likelihood_profile.opju` project and apply a filter on the error column. Set the filter to <0.2 initially
23. Chek 1) whether the filtered region size is comparable to the total explored region, 2) there are enough points populating the plot
24. Decide whether to expand or shrink the explored region and whether more points are needed to fill out the interest region
25. Iterate over steps 15-22 until you get at least two filteres regions, for example, E<0.1 and E<0.2 with enough points to showcase that the region shrinks as the error acceptance decreases  
``

#### Local Sensitivity Analysis

1.  Access a Slurm HPC cluster
2.  Transfer the cluster scripts in `cluster -> LSA`  to your cluster folder in the HPC
3.  You may need to set up Anaconda/Miniconda in your HPC cluster
4.  Install CompuCell3D in your cluster using the conda installation (see previous section or follow the steps in the [Official Website](https://compucell3d.org/SrcBin))
5.  Activate your CC3D environment
6.  Adapt the scripts to your system, such as paths, user name etc:
    ```bash
    # Open the 'script.py' file

    # Set the wanted number of replicas, total number of parameter sets, or leave as is
    replicas = 200
    N = 1000

    # Set the maximum number of jobs. This depends on the specifics of your cluster
    MAX_JOBS = 9000

    # Set the multiplier, or leave as is. This value sets by how much each parameter is going to change for comparison to default
    multiplier = 1.5
    ```
7. In the generate_sbatch_string function, change output, email and path to script
8. Run the script by calling `python script.py
9. Collect the output files and run the `LSA.py` script on the output files
11. Open the `LSA_final.opju` project (OriginLab required) from the `Results_and_Plots` folder
12. Transfer the averaged table to the `LSA_final.opju` project and plot the metrics versus the parameter names using bar plots
13. Alternatively, before plotting the metrics, you can order the output metric from low to high, and use the same ordering for the other plots, fixing the parameter names in place and making it easier to compare between different cases 

#### Exploring scenarios

There are a couple of scenarios you can explore with the simulation. You can do this with the online version too
* **NS1 deletion**: in the parameter table, you can emulate an NS1 deletion scenario by setting the defect probability of NS1 to 1. In the final plots, you will notice lower IFN production compared to the default parameter case
* **IFN pretreatment / high baseline IFN secretion**: in the parameter table, set the baseline secretion of IFN to a higher value such as 1
* **MOI regimes**: try a very high MOI such as MOI=1. You will see faster peaks of virus and IFN production, and the IFN production at high MOI with NS1 competent viruses is much lower than for lower MOI

---

## Citation

> Waiting for review

If you use this model or the associated code in your research, please cite both the manuscript and the Zenodo software archive.

### Manuscript

```bibtex
@article{handle,
  author = {Authors},
  title = {{Title}},
  journal = {journal},
  year = {2026},
  doi = {0000},
  url = {[website)}
}
```

### Software

```bibtex
@software{handle,
  author = {Authors},
  title = {{Title}},
  month = June,
  year = 2026,
  publisher = {repository},
  version = {version},
  doi = {0000},
  url = {[website)}
}
