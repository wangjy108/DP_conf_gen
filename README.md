# sample usage for base util

**requiest**

pre-requst: rdkit logging scipy

```bash

git clone https://github.com/wangjy108/base_util.git ${INSTALL_DIR}

cd ${INSTALL_DIR}                   

python setup.py install    # run

```

**Docker**
```bash
docker pull registry.dp.tech/dptech/prod-1364/strain:run0.0.2
docker run -dit -v ${PWD}:/data --name strain registry.dp.tech/dptech/prod-1364/strain:run0.0.2 /bin/bash  # set image
docker attach strain
cd /data
# initiate config.in file
python ${each}                                              # run

```