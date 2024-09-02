#!/bin/bash

num_cores=4
res_sp=5

# Definir el rango de núcleos a utilizar (por ejemplo, de 1 a 8 núcleos)
for res_tmp in $(seq 1 1 4)
do
    # Ejecutar la simulación con la cantidad actual de núcleos y capturar la salida en una variable
    output=$(mpirun -np $num_cores ./cavity3d $res_sp $res_tmp)

    # Filtrar la línea que contiene "realTime" y "cpuTime"
    filtered_output=$(echo "$output" | grep "realTime=")

    # Extraer valores usando awk
    real_time=$(echo "$filtered_output" | awk -F'[=;]' '{print $2}')
    cpu_time=$(echo "$filtered_output" | awk -F'[=;]' '{print $4}')

    # Imprimir los valores en la consola junto con el número de núcleos
    echo "M: $res_tmp"
    echo "realTime=${real_time}"
    echo "cpuTime=${cpu_time}"
    echo

    # Guardar los valores en un archivo junto con el número de núcleos
    echo "$res_tmp,${real_time},${cpu_time}" >> resultados.txt
done
