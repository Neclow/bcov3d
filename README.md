# Recombination of binding domains shapes structural evolution of Betacoronavirus spike proteins

## Installation

Note: this project was implemented using Linux 64-bit.

### Dependencies

This project uses [pixi](https://pixi.sh/) for dependency management.

To install the dependencies, run:

```bash
# Install conda/pypi dependencies
pixi install
# Set up & build other dependencies (envalysis, FAMSA3di)
pixi run post_install
```

#### InterProScan

This project uses Docker to run InterProScan (v5.76-107.0), for which you will need to have superuser access. Installing the InterProScan database requires 35G of space.

Follow [these instructions](https://docs.docker.com/engine/install/ubuntu) to install Docker. See [InterProScan usage instructions](https://interproscan-docs.readthedocs.io/en/v5/HowToUseViaContainer.html) for more details.

#### Famsa3di

We use a custom fork of [FAMSA](https://github.com/refresh-bio/FAMSA) to make it compatible with 3di, drawing from the [3diphy](https://github.com/nmatzke/3diphy) project.

To add the famsa3di submodule manually, use:

```bash
git submodule add https://github.com/Neclow/famsa3di extern/famsa3di
cd extern/famsa3di
git submodule update --init --recursive
```

### Data

TBD

