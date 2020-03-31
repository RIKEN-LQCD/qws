#!/bin/bash -x

#for j in utofu utofu_threaded mpi_rankmap;do
for j in utofu ;do
for i in all_calc overlapped send send_post; do
  make clean
  make profiler=timing2 target=$i rdma=$j >& make.$i.$j.log
  mv main main_${i}_${j}
done
done
