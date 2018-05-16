rm -rf roofit_input_C.so
rm -rf roofit_input_C.d
rm -rf result.log

root -l -b  < x_roofit_input.C  &> result.log &

rm -rf roofit_input_C.so
rm -rf roofit_input_C.d
