Build and locally test the website
==================================

You need to have docker installed.

### Step 1: Build

To build

	docker build -t bioconductor_website:latest

### Step 2: Run 

To test if your changes have made it into

	docker run -p 3000:3000 bioc_website:latest

### Step 3: Open your browser http://localhost:3000

To see the latest website, open http://localhost:3000


