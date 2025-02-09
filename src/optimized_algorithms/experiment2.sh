#!/bin/bash

# echo "Hello, GeeksforGeeks"
# for ((domain_size_x=256;domain_size_x<=512;domain_size_x=domain_size_x*2))
#     do
#     for ((domain_size_y=256;domain_size_y<=512;domain_size_y=domain_size_y*2))
#         do
#         for ((domain_size_z=256;domain_size_z<=512;domain_size_z=domain_size_z*2))
#             do
#             for ((nbins=2;nbins<=6000;nbins=nbins*2))
#                 do

#                 # nbins=2
#                 # domain_size=256
#                 max_grid_size=$((domain_size_x/nbins))
#                 grid_size=$((max_grid_size/2))
                
#                 # boxes=8
#                 for ((i=max_grid_size; i>=grid_size && grid_size > 0; i=i/2))
#                 do
#                     for ((j=i;  j>=grid_size; j=j/2))
#                     do
#                         for((k=j; k>=grid_size; k=k/2))
#                         do
#                             echo "i=$i,j=$j,k=$k"
#                             boxes=$(($((domain_size/i))*$((domain_size/j))*$((domain_size/k))))
#                             echo "nbins=$nbins"
#                             echo "boxes=$boxes"
#                             echo "max_grid_size=$max_grid_size"
#                             echo  "grid_size=$grid_size"
#                             cmd='./main3d.gnu.x86-milan.TPROF.ex inputs  nbins='$nbins' domain="('$domain_size_x,$domain_size_y,$domain_size_z')" max_grid_size="('$i,$j,$k')" >output/'$domain_size'_'$nbins'_'$boxes'_output.txt'
#                             echo $cmd
#                             # eval $cmd
#                         done
                        
#                     done
#                 done
#             done
#         done
#     done
# done
# ./main3d.gnu.x86-milan.TPROF.OMP.ex inputs 
# ./main3d.gnu.x86-milan.TPROF.ex inputs 
# ./main3d.gnu.x86-milan.TPROF.ex inputs  nnodes=128 max_grid_size="(32,16,16)" >output_withfixedseed/128_2048_output.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs  nbins=2 max_grid_size="(128,128,64)" >output/2_16_combined_output_output.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs  nbins=2 max_grid_size="(128,64,64)" >output/2_32_output.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs  nbins=2 max_grid_size="(64,64,64)" >output/2_64_combined_output.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs  nbins=2 max_grid_size="(64,64,32)" >output/2_128_combined_output.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs  nbins=2 max_grid_size="(64,32,32)" >output/2_256_combined_output.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs  nbins=2 max_grid_size="(32,32,32)" >output/2_512_combined_output.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs  nbins=2 max_grid_size="(32,32,16)" >output/2_1024_output.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs  nbins=2 max_grid_size="(32,16,16)" >output/2_2048_output.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs  nbins=2 max_grid_size="(16,16,16)" >output/2_4096_output.txt



#Uncomment the test case which you want to run

# cd src

# 1 node 4 ranks 4 boxes per rank


########## Best Case ############

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=1 domain="(256,256,256)" max_grid_size="(128,128,64)" >../../output/run_250/4_4_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=1 domain="(256,256,256)" max_grid_size="(128,64,64)" >../../output/run_250/4_8_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=1 domain="(256,256,256)" max_grid_size="(64,64,64)" >../../output/run_250/4_16_output_best.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2 domain="(256,256,256)" max_grid_size="(128,64,64)" >../../output/run_250/8_4_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2 domain="(256,256,256)" max_grid_size="(64,64,64)" >../../output/run_250/8_8_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2 domain="(256,256,256)" max_grid_size="(64,64,32)" >../../output/run_250/8_16_output_best.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=4 domain="(256,256,256)" max_grid_size="(64,64,64)" >../../output/run_250/16_4_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=4 domain="(256,256,256)" max_grid_size="(64,64,32)" >../../output/run_250/16_8_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=4 domain="(256,256,256)" max_grid_size="(64,32,32)" >../../output/run_250/16_16_output_best.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=8 domain="(256,256,256)" max_grid_size="(64,64,32)" >../../output/run_250/32_4_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=8 domain="(256,256,256)" max_grid_size="(64,32,32)" >../../output/run_250/32_8_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=8 domain="(256,256,256)" max_grid_size="(32,32,32)" >../../output/run_250/32_16_output_best.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=16 domain="(256,256,256)" max_grid_size="(64,32,32)" >../../output/run_250/64_4_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=16 domain="(256,256,256)" max_grid_size="(32,32,32)" >../../output/run_250/64_8_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=16 domain="(256,256,256)" max_grid_size="(32,32,16)" >../../output/run_250/64_16_output_best.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=32 domain="(256,256,256)" max_grid_size="(32,32,32)" >../../output/run_250/128_4_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=32 domain="(256,256,256)" max_grid_size="(32,32,16)" >../../output/run_250/128_8_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=32 domain="(256,256,256)" max_grid_size="(32,16,16)" >../../output/run_250/128_16_output_best.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=64 domain="(256,256,256)" max_grid_size="(32,32,16)" >../../output/run_250/256_4_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=64 domain="(256,256,256)" max_grid_size="(32,16,16)" >../../output/run_250/256_8_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=64 domain="(256,256,256)" max_grid_size="(16,16,16)" >../../output/run_250/256_16_output_best.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=128 domain="(256,256,256)" max_grid_size="(32,16,16)" >../../output/run_250/512_4_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=128 domain="(256,256,256)" max_grid_size="(16,16,16)" >../../output/run_250/512_8_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=128 domain="(256,256,256)" max_grid_size="(16,16,8)" >../../output/run_250/512_16_output_best.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=256 domain="(256,256,256)" max_grid_size="(16,16,16)" >../../output/run_250/1024_4_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=256 domain="(256,256,256)" max_grid_size="(16,16,8)" >../../output/run_250/1024_8_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=256 domain="(256,256,256)" max_grid_size="(16,8,8)" >../../output/run_250/1024_16_output_best.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=512 domain="(256,256,256)" max_grid_size="(16,16,8)" >../../output/run_250/2048_4_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=512 domain="(256,256,256)" max_grid_size="(16,8,8)" >../../output/run_250/2048_8_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=512 domain="(256,256,256)" max_grid_size="(8,8,8)" >../../output/run_250/2048_16_output_best.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=1024 domain="(256,256,256)" max_grid_size="(16,8,8)" >../../output/run_250/4096_4_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=1024 domain="(256,256,256)" max_grid_size="(8,8,8)" >../../output/run_250/4096_8_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=1024 domain="(256,256,256)" max_grid_size="(8,8,4)" >../../output/run_250/4096_16_output_best.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2048 domain="(256,256,256)" max_grid_size="(8,8,8)" >../../output/run_250/8192_4_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2048 domain="(256,256,256)" max_grid_size="(8,8,4)" >../../output/run_250/8192_8_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2048 domain="(256,256,256)" max_grid_size="(8,4,4)" >../../output/run_250/8192_16_output_best.txt




########## Avg Case ############



# ./main3d.gnu.x86-milan.TPROF.ex inputs1 nnodes=1 domain="(256,256,256)" max_grid_size="(128,128,64)" >../../output/run_250/4_4_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=1 domain="(256,256,256)" max_grid_size="(128,64,64)" >../../output/run_250/4_8_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=1 domain="(256,256,256)" max_grid_size="(64,64,64)" >../../output/run_250/4_16_output_avg.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2 domain="(256,256,256)" max_grid_size="(128,64,64)" >../../output/run_250/8_4_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2 domain="(256,256,256)" max_grid_size="(64,64,64)" >../../output/run_250/8_8_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2 domain="(256,256,256)" max_grid_size="(64,64,32)" >../../output/run_250/8_16_output_avg.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=4 domain="(256,256,256)" max_grid_size="(64,64,64)" >../../output/run_250/16_4_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=4 domain="(256,256,256)" max_grid_size="(64,64,32)" >../../output/run_250/16_8_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=4 domain="(256,256,256)" max_grid_size="(64,32,32)" >../../output/run_250/16_16_output_avg.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=8 domain="(256,256,256)" max_grid_size="(64,64,32)" >../../output/run_250/32_4_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=8 domain="(256,256,256)" max_grid_size="(64,32,32)" >../../output/run_250/32_8_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=8 domain="(256,256,256)" max_grid_size="(32,32,32)" >../../output/run_250/32_16_output_avg.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=16 domain="(256,256,256)" max_grid_size="(64,32,32)" >../../output/run_250/64_4_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=16 domain="(256,256,256)" max_grid_size="(32,32,32)" >../../output/run_250/64_8_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=16 domain="(256,256,256)" max_grid_size="(32,32,16)" >../../output/run_250/64_16_output_avg.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=32 domain="(256,256,256)" max_grid_size="(32,32,32)" >../../output/run_250/128_4_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=32 domain="(256,256,256)" max_grid_size="(32,32,16)" >../../output/run_250/128_8_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=32 domain="(256,256,256)" max_grid_size="(32,16,16)" >../../output/run_250/128_16_output_avg.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=64 domain="(256,256,256)" max_grid_size="(32,32,16)" >../../output/run_250/256_4_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=64 domain="(256,256,256)" max_grid_size="(32,16,16)" >../../output/run_250/256_8_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=64 domain="(256,256,256)" max_grid_size="(16,16,16)" >../../output/run_250/256_16_output_avg.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=128 domain="(256,256,256)" max_grid_size="(32,16,16)" >../../output/run_250/512_4_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=128 domain="(256,256,256)" max_grid_size="(16,16,16)" >../../output/run_250/512_8_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=128 domain="(256,256,256)" max_grid_size="(16,16,8)" >../../output/run_250/512_16_output_avg.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=256 domain="(256,256,256)" max_grid_size="(16,16,16)" >../../output/run_250/1024_4_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=256 domain="(256,256,256)" max_grid_size="(16,16,8)" >../../output/run_250/1024_8_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=256 domain="(256,256,256)" max_grid_size="(16,8,8)" >../../output/run_250/1024_16_output_avg.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=512 domain="(256,256,256)" max_grid_size="(16,16,8)" >../../output/run_250/2048_4_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=512 domain="(256,256,256)" max_grid_size="(16,8,8)" >../../output/run_250/2048_8_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=512 domain="(256,256,256)" max_grid_size="(8,8,8)" >../../output/run_250/2048_16_output_avg.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=1024 domain="(256,256,256)" max_grid_size="(16,8,8)" >../../output/run_250/4096_4_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=1024 domain="(256,256,256)" max_grid_size="(8,8,8)" >../../output/run_250/4096_8_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=1024 domain="(256,256,256)" max_grid_size="(8,8,4)" >../../output/run_250/4096_16_output_avg.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2048 domain="(256,256,256)" max_grid_size="(8,8,8)" >../../output/run_250/8192_4_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2048 domain="(256,256,256)" max_grid_size="(8,8,4)" >../../output/run_250/8192_8_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2048 domain="(256,256,256)" max_grid_size="(8,4,4)" >../../output/run_250/8192_16_output_avg.txt



########## Worst Case ############



# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=1 domain="(256,256,256)" max_grid_size="(128,128,64)" >../../output/run_250/4_4_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=1 domain="(256,256,256)" max_grid_size="(128,64,64)" >../../output/run_250/4_8_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=1 domain="(256,256,256)" max_grid_size="(64,64,64)" >../../output/run_250/4_16_output_worst.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=2 domain="(256,256,256)" max_grid_size="(128,64,64)" >../../output/run_250/8_4_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=2 domain="(256,256,256)" max_grid_size="(64,64,64)" >../../output/run_250/8_8_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=2 domain="(256,256,256)" max_grid_size="(64,64,32)" >../../output/run_250/8_16_output_worst.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=4 domain="(256,256,256)" max_grid_size="(64,64,64)" >../../output/run_250/16_4_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=4 domain="(256,256,256)" max_grid_size="(64,64,32)" >../../output/run_250/16_8_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=4 domain="(256,256,256)" max_grid_size="(64,32,32)" >../../output/run_250/16_16_output_worst.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=8 domain="(256,256,256)" max_grid_size="(64,64,32)" >../../output/run_250/32_4_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=8 domain="(256,256,256)" max_grid_size="(64,32,32)" >../../output/run_250/32_8_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=8 domain="(256,256,256)" max_grid_size="(32,32,32)" >../../output/run_250/32_16_output_worst.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=16 domain="(256,256,256)" max_grid_size="(64,32,32)" >../../output/run_250/64_4_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=16 domain="(256,256,256)" max_grid_size="(32,32,32)" >../../output/run_250/64_8_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=16 domain="(256,256,256)" max_grid_size="(32,32,16)" >../../output/run_250/64_16_output_worst.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=32 domain="(256,256,256)" max_grid_size="(32,32,32)" >../../output/run_250/128_4_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=32 domain="(256,256,256)" max_grid_size="(32,32,16)" >../../output/run_250/128_8_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=32 domain="(256,256,256)" max_grid_size="(32,16,16)" >../../output/run_250/128_16_output_worst.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=64 domain="(256,256,256)" max_grid_size="(32,32,16)" >../../output/run_250/256_4_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=64 domain="(256,256,256)" max_grid_size="(32,16,16)" >../../output/run_250/256_8_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=64 domain="(256,256,256)" max_grid_size="(16,16,16)" >../../output/run_250/256_16_output_worst.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=128 domain="(256,256,256)" max_grid_size="(32,16,16)" >../../output/run_250/512_4_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=128 domain="(256,256,256)" max_grid_size="(16,16,16)" >../../output/run_250/512_8_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=128 domain="(256,256,256)" max_grid_size="(16,16,8)" >../../output/run_250/512_16_output_worst.txt

./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=256 domain="(256,256,256)" max_grid_size="(16,16,16)" >../../output/run_250/1024_4_output_worst.txt
./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=256 domain="(256,256,256)" max_grid_size="(16,16,8)" >../../output/run_250/1024_8_output_worst.txt
./main3d.gnu.x86-milan.TPROF.ex inputs2 nnodes=256 domain="(256,256,256)" max_grid_size="(16,8,8)" >../../output/run_250/1024_16_output_worst.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=512 domain="(256,256,256)" max_grid_size="(16,16,8)" >../../output/run_250/2048_4_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=512 domain="(256,256,256)" max_grid_size="(16,8,8)" >../../output/run_250/2048_8_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=512 domain="(256,256,256)" max_grid_size="(8,8,8)" >../../output/run_250/2048_16_output_worst.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=1024 domain="(256,256,256)" max_grid_size="(16,8,8)" >../../output/run_250/4096_4_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=1024 domain="(256,256,256)" max_grid_size="(8,8,8)" >../../output/run_250/4096_8_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=1024 domain="(256,256,256)" max_grid_size="(8,8,4)" >../../output/run_250/4096_16_output_worst.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2048 domain="(256,256,256)" max_grid_size="(8,8,8)" >../../output/run_250/8192_4_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2048 domain="(256,256,256)" max_grid_size="(8,8,4)" >../../output/run_250/8192_8_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2048 domain="(256,256,256)" max_grid_size="(8,4,4)" >../../output/run_250/8192_16_output_worst.txt






















######### old ones #############

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2 domain="(256,256,256)" max_grid_size="(128,128,128)" >../../output/run_250/2_4_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2 domain="(256,256,256)" max_grid_size="(128,128,64)" >../../output/run_250/2_8_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2 domain="(256,256,256)" max_grid_size="(128,64,64)" >../../output/run_250/2_16_output_best.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=4 domain="(256,256,256)" max_grid_size="(128,128,64)" >../../output/run_250/4_4_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=4 domain="(256,256,256)" max_grid_size="(128,64,64)" >../../output/run_250/4_8_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=4 domain="(256,256,256)" max_grid_size="(64,64,64)" >../../output/run_250/4_16_output_best.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=8 domain="(256,256,256)" max_grid_size="(128,64,64)" >../../output/run_250/8_4_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=8 domain="(256,256,256)" max_grid_size="(64,64,64)" >../../output/run_250/8_8_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=8 domain="(256,256,256)" max_grid_size="(64,64,32)" >../../output/run_250/8_16_output_best.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=16 domain="(256,256,256)" max_grid_size="(64,64,64)" >../../output/run_250/16_4_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=16 domain="(256,256,256)" max_grid_size="(64,64,32)" >../../output/run_250/16_8_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=16 domain="(256,256,256)" max_grid_size="(64,32,32)" >../../output/run_250/16_16_output_best.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=32 domain="(256,256,256)" max_grid_size="(64,64,32)" >../../output/run_250/32_4_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=32 domain="(256,256,256)" max_grid_size="(64,32,32)" >../../output/run_250/32_8_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=32 domain="(256,256,256)" max_grid_size="(32,32,32)" >../../output/run_250/32_16_output_best.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=64 domain="(256,256,256)" max_grid_size="(64,32,32)" >../../output/run_250/64_4_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=64 domain="(256,256,256)" max_grid_size="(32,32,32)" >../../output/run_250/64_8_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=64 domain="(256,256,256)" max_grid_size="(32,32,16)" >../../output/run_250/64_16_output_best.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=128 domain="(256,256,256)" max_grid_size="(32,32,32)" >../../output/run_250/128_4_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=128 domain="(256,256,256)" max_grid_size="(32,32,16)" >../../output/run_250/128_8_output_best.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=128 domain="(256,256,256)" max_grid_size="(32,16,16)" >../../output/run_250/128_16_output_best.txt



# #########################

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2 domain="(256,256,256)" max_grid_size="(128,128,128)" >../../output/run_250/2_4_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2 domain="(256,256,256)" max_grid_size="(128,128,64)" >../../output/run_250/2_8_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2 domain="(256,256,256)" max_grid_size="(128,64,64)" >../../output/run_250/2_16_output_avg.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=4 domain="(256,256,256)" max_grid_size="(128,128,64)" >../../output/run_250/4_4_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=4 domain="(256,256,256)" max_grid_size="(128,64,64)" >../../output/run_250/4_8_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=4 domain="(256,256,256)" max_grid_size="(64,64,64)" >../../output/run_250/4_16_output_avg.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=8 domain="(256,256,256)" max_grid_size="(128,64,64)" >../../output/run_250/8_4_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=8 domain="(256,256,256)" max_grid_size="(64,64,64)" >../../output/run_250/8_8_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=8 domain="(256,256,256)" max_grid_size="(64,64,32)" >../../output/run_250/8_16_output_avg.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=16 domain="(256,256,256)" max_grid_size="(64,64,64)" >../../output/run_250/16_4_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=16 domain="(256,256,256)" max_grid_size="(64,64,32)" >../../output/run_250/16_8_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=16 domain="(256,256,256)" max_grid_size="(64,32,32)" >../../output/run_250/16_16_output_avg.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=32 domain="(256,256,256)" max_grid_size="(64,64,32)" >../../output/run_250/32_4_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=32 domain="(256,256,256)" max_grid_size="(64,32,32)" >../../output/run_250/32_8_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=32 domain="(256,256,256)" max_grid_size="(32,32,32)" >../../output/run_250/32_16_output_avg.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=64 domain="(256,256,256)" max_grid_size="(64,32,32)" >../../output/run_250/64_4_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=64 domain="(256,256,256)" max_grid_size="(32,32,32)" >../../output/run_250/64_8_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=64 domain="(256,256,256)" max_grid_size="(32,32,16)" >../../output/run_250/64_16_output_avg.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=128 domain="(256,256,256)" max_grid_size="(32,32,32)" >../../output/run_250/128_4_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=128 domain="(256,256,256)" max_grid_size="(32,32,16)" >../../output/run_250/128_8_output_avg.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=128 domain="(256,256,256)" max_grid_size="(32,16,16)" >../../output/run_250/128_16_output_avg.txt

# # #########################

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2 domain="(256,256,256)" max_grid_size="(128,128,128)" >../../output/run_250/2_4_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2 domain="(256,256,256)" max_grid_size="(128,128,64)" >../../output/run_250/2_8_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=2 domain="(256,256,256)" max_grid_size="(128,64,64)" >../../output/run_250/2_16_output_worst.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=4 domain="(256,256,256)" max_grid_size="(128,128,64)" >../../output/run_250/4_4_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=4 domain="(256,256,256)" max_grid_size="(128,64,64)" >../../output/run_250/4_8_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=4 domain="(256,256,256)" max_grid_size="(64,64,64)" >../../output/run_250/4_16_output_worst.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=8 domain="(256,256,256)" max_grid_size="(128,64,64)" >../../output/run_250/8_4_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=8 domain="(256,256,256)" max_grid_size="(64,64,64)" >../../output/run_250/8_8_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=8 domain="(256,256,256)" max_grid_size="(64,64,32)" >../../output/run_250/8_16_output_worst.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=16 domain="(256,256,256)" max_grid_size="(64,64,64)" >../../output/run_250/16_4_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=16 domain="(256,256,256)" max_grid_size="(64,64,32)" >../../output/run_250/16_8_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=16 domain="(256,256,256)" max_grid_size="(64,32,32)" >../../output/run_250/16_16_output_worst.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=32 domain="(256,256,256)" max_grid_size="(64,64,32)" >../../output/run_250/32_4_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=32 domain="(256,256,256)" max_grid_size="(64,32,32)" >../../output/run_250/32_8_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=32 domain="(256,256,256)" max_grid_size="(32,32,32)" >../../output/run_250/32_16_output_worst.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=64 domain="(256,256,256)" max_grid_size="(64,32,32)" >../../output/run_250/64_4_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=64 domain="(256,256,256)" max_grid_size="(32,32,32)" >../../output/run_250/64_8_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=64 domain="(256,256,256)" max_grid_size="(32,32,16)" >../../output/run_250/64_16_output_worst.txt

# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=128 domain="(256,256,256)" max_grid_size="(32,32,32)" >../../output/run_250/128_4_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=128 domain="(256,256,256)" max_grid_size="(32,32,16)" >../../output/run_250/128_8_output_worst.txt
# ./main3d.gnu.x86-milan.TPROF.ex inputs nnodes=128 domain="(256,256,256)" max_grid_size="(32,16,16)" >../../output/run_250/128_16_output_worst.txt


