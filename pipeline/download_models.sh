MODEL_DIR="data/models"
mkdir -p ${MODEL_DIR}

wget -O ${MODEL_DIR}/Q.3DI.AF https://edmond.mpg.de/api/access/datafile/311466
wget -O ${MODEL_DIR}/Q.3DI.LLM https://edmond.mpg.de/api/access/datafile/311467
wget -O ${MODEL_DIR}/VK23 https://raw.githubusercontent.com/nmatzke/3diphy/c6e1e6cc334776b69dd20f9a0a5a306bf7cb98c5/3DI_substmat/3DI
wget -qO- https://raw.githubusercontent.com/steineggerlab/foldseek/refs/tags/10-941cd33/data/mat3di.out | tail -n +4 > ${MODEL_DIR}/VK23.scoring

# IQTREE3 has a bug where it fails to read substitution matrices with paths, so we also create symlinks in the current directory
ln -sf ${MODEL_DIR}/Q.3DI.AF Q.3DI.AF
ln -sf ${MODEL_DIR}/Q.3DI.LLM Q.3DI.LLM
ln -sf ${MODEL_DIR}/VK23 VK23
ln -sf ${MODEL_DIR}/VK23.scoring VK23.scoring