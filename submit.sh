for i in $(ls ./data/xyz/*.xyz); do
    tmp=${i%.*}
    name=${tmp:11}
    echo $name

    rm -rf ./work/$name/*
    mkdir ./work/$name
    
    cp ./data/xyz/$name.xyz ./work/$name/py-ch3-cl.xyz
    cp ./src/run.sh  ./work/$name/run.sh
    cp ./src/main.py ./work/$name/main.py

    # cd ./work/$name;
    # sbatch run.sh;
    # cd -;
done