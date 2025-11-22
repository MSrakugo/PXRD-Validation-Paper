# Programs and Datasets for Validating PXRD mineral modes

This repository contains programs and datasets for validating PXRD mineral modes. The preprint or published paper URL will be available soon.

## Getting Started

To start this project, you need to set up and activate the Python **3.11** environment using **Poetry**, and then install the necessary dependencies.

### **Step 1: Specify the Python Environment**

Use the `poetry env use` command to tell Poetry which Python version to use for the project.

```bash
poetry env use 3.11
```

### **Step 2: Activate the Environment**

Activate the environment you've just configured. This ensures that the project's dependencies are run within an isolated virtual environment.

```bash
poetry env activate
```

### **Step 3: Install Dependencies**

Once the environment is active, run the following command to install all the required project dependencies as listed in `pyproject.toml`.

```bash
poetry install
```

### **Next step: Make figures with notebooks**

Currently, you can run programs written in notebooks. With a dataset in the folder named "0_Data," you can run the following notebooks: `2_Program/0_2_Visualize_Data_251110.ipynb`. If you cannot run a notebook, you can view all notebook outputs.

## Author

Satoshi Matsuno, Tohoku university, Japan