#./analysisbash,sh


echo "Enter run type: SRC or MF"
read run

echo "Enter pass: pass1, pass2, pass3, or pass4"
read pass

echo "Enter coin/sing: coin or sing"
read sorc

file="omniSkims.C"
location="pwd"


pathingBe9="cd Src/pass4/Be9/Coin/"
pathingB10="cd ../../B10/Coin/"
pathingB11="cd ../../B11/Coin/"
pathingC12="cd ../../C12/Coin/"
pathingCa40="cd ../../Ca40/Coin/"
pathingCa48="cd ../../Ca48/Coin/"
pathingFe54="cd ../../Fe54/Coin/"
pathingAu197="cd ../../Au197/Coin/"
pathingD2="cd ../../D2/Coin/"

if [ $run == 'SRC' ]
	then
		if [ $pass == 'pass1' ]
			then
				if [ $sorc == 'coin' ]
					then
						pathingBe9="cd Src/pass1/Be9/Coin/"
						pathingB10="cd ../../B10/Coin/"
						pathingB11="cd ../../B11/Coin/"
						pathingC12="cd ../../C12/Coin/"
						pathingCa40="cd ../../Ca40/Coin/"
						pathingCa48="cd ../../Ca48/Coin/"
						pathingFe54="cd ../../Fe54/Coin/"
						pathingAu197="cd ../../Au197/Coin/"
						pathingD2="cd ../../D2/Coin/"
				elif [ $sorc == 'sing' ]
					then
						pathingBe9="cd Src/pass1/Be9/Sing/"
						pathingB10="cd ../../B10/Sing/"
						pathingB11="cd ../../B11/Sing/"
						pathingC12="cd ../../C12/Sing/"
						pathingCa40="cd ../../Ca40/Sing/"
						pathingCa48="cd ../../Ca48/Sing/"
						pathingFe54="cd ../../Fe54/Sing/"
						pathingAu197="cd ../../Au197/Sing/"
						pathingD2="cd ../../D2/Sing/"
				fi
		elif [ $pass == 'pass2' ]
			then
				if [ $sorc == 'coin' ]
					then
						pathingBe9="cd Src/pass2/Be9/Coin/"
						pathingB10="cd ../../B10/Coin/"
						pathingB11="cd ../../B11/Coin/"
						pathingC12="cd ../../C12/Coin/"
						pathingCa40="cd ../../Ca40/Coin/"
						pathingCa48="cd ../../Ca48/Coin/"
						pathingFe54="cd ../../Fe54/Coin/"
						pathingAu197="cd ../../Au197/Coin/"
						pathingD2="cd ../../D2/Coin/"
				elif [ $sorc == 'sing' ]
					then
						pathingBe9="cd Src/pass2/Be9/Sing/"
						pathingB10="cd ../../B10/Sing/"
						pathingB11="cd ../../B11/Sing/"
						pathingC12="cd ../../C12/Sing/"
						pathingCa40="cd ../../Ca40/Sing/"
						pathingCa48="cd ../../Ca48/Sing/"
						pathingFe54="cd ../../Fe54/Sing/"
						pathingAu197="cd ../../Au197/Sing/"
						pathingD2="cd ../../D2/Sing/"
				fi
		elif [ $pass == 'pass3' ]
			then
				if [ $sorc == 'coin' ]
					then
						pathingBe9="cd Src/pass3/Be9/Coin/"
						pathingB10="cd ../../B10/Coin/"
						pathingB11="cd ../../B11/Coin/"
						pathingC12="cd ../../C12/Coin/"
						pathingCa40="cd ../../Ca40/Coin/"
						pathingCa48="cd ../../Ca48/Coin/"
						pathingFe54="cd ../../Fe54/Coin/"
						pathingAu197="cd ../../Au197/Coin/"
						pathingD2="cd ../../D2/Coin/"
				elif [ $sorc == 'sing' ]
					then
						pathingBe9="cd Src/pass3/Be9/Sing/"
						pathingB10="cd ../../B10/Sing/"
						pathingB11="cd ../../B11/Sing/"
						pathingC12="cd ../../C12/Sing/"
						pathingCa40="cd ../../Ca40/Sing/"
						pathingCa48="cd ../../Ca48/Sing/"
						pathingFe54="cd ../../Fe54/Sing/"
						pathingAu197="cd ../../Au197/Sing/"
						pathingD2="cd ../../D2/Sing/"
				fi
		elif [ $pass == 'pass4' ]
			then
				if [ $sorc == 'coin' ]
					then
						pathingBe9="cd Src/pass4/Be9/Coin/"
						pathingB10="cd ../../B10/Coin/"
						pathingB11="cd ../../B11/Coin/"
						pathingC12="cd ../../C12/Coin/"
						pathingCa40="cd ../../Ca40/Coin/"
						pathingCa48="cd ../../Ca48/Coin/"
						pathingFe54="cd ../../Fe54/Coin/"
						pathingAu197="cd ../../Au197/Coin/"
						pathingD2="cd ../../D2/Coin/"
				elif [ $sorc == 'sing' ]
					then
						pathingBe9="cd Src/pass4/Be9/Sing/"
						pathingB10="cd ../../B10/Sing/"
						pathingB11="cd ../../B11/Sing/"
						pathingC12="cd ../../C12/Sing/"
						pathingCa40="cd ../../Ca40/Sing/"
						pathingCa48="cd ../../Ca48/Sing/"
						pathingFe54="cd ../../Fe54/Sing/"
						pathingAu197="cd ../../Au197/Sing/"
						pathingD2="cd ../../D2/Sing/"
				fi
		fi
elif [ $run == 'MF' ]
	then
		if [ $pass == 'pass1' ]
			then
				if [ $sorc == 'coin' ]
					then
						pathingBe9="cd Mf/pass1/Be9/Coin/"
						pathingB10="cd ../../B10/Coin/"
						pathingB11="cd ../../B11/Coin/"
						pathingC12="cd ../../C12/Coin/"
						pathingCa40="cd ../../Ca40/Coin/"
						pathingCa48="cd ../../Ca48/Coin/"
						pathingFe54="cd ../../Fe54/Coin/"
						pathingAu197="cd ../../Au197/Coin/"
						pathingD2="cd ../../D2/Coin/"
				elif [ $sorc == 'sing' ]
					then
						pathingBe9="cd Mf/pass1/Be9/Sing/"
						pathingB10="cd ../../B10/Sing/"
						pathingB11="cd ../../B11/Sing/"
						pathingC12="cd ../../C12/Sing/"
						pathingCa40="cd ../../Ca40/Sing/"
						pathingCa48="cd ../../Ca48/Sing/"
						pathingFe54="cd ../../Fe54/Sing/"
						pathingAu197="cd ../../Au197/Sing/"
						pathingD2="cd ../../D2/Sing/"
				fi
		elif [ $pass == 'pass2' ]
			then
				if [ $sorc == 'coin' ]
					then
						pathingBe9="cd Mf/pass2/Be9/Coin/"
						pathingB10="cd ../../B10/Coin/"
						pathingB11="cd ../../B11/Coin/"
						pathingC12="cd ../../C12/Coin/"
						pathingCa40="cd ../../Ca40/Coin/"
						pathingCa48="cd ../../Ca48/Coin/"
						pathingFe54="cd ../../Fe54/Coin/"
						pathingAu197="cd ../../Au197/Coin/"
						pathingD2="cd ../../D2/Coin/"
				elif [ $sorc == 'sing' ]
					then
						pathingBe9="cd Mf/pass2/Be9/Sing/"
						pathingB10="cd ../../B10/Sing/"
						pathingB11="cd ../../B11/Sing/"
						pathingC12="cd ../../C12/Sing/"
						pathingCa40="cd ../../Ca40/Sing/"
						pathingCa48="cd ../../Ca48/Sing/"
						pathingFe54="cd ../../Fe54/Sing/"
						pathingAu197="cd ../../Au197/Sing/"
						pathingD2="cd ../../D2/Sing/"
				fi
		elif [ $pass == 'pass3' ]
			then
				if [ $sorc == 'coin' ]
					then
						pathingBe9="cd Mf/pass3/Be9/Coin/"
						pathingB10="cd ../../B10/Coin/"
						pathingB11="cd ../../B11/Coin/"
						pathingC12="cd ../../C12/Coin/"
						pathingCa40="cd ../../Ca40/Coin/"
						pathingCa48="cd ../../Ca48/Coin/"
						pathingFe54="cd ../../Fe54/Coin/"
						pathingAu197="cd ../../Au197/Coin/"
						pathingD2="cd ../../D2/Coin/"
				elif [ $sorc == 'sing' ]
					then
						pathingBe9="cd Mf/pass3/Be9/Sing/"
						pathingB10="cd ../../B10/Sing/"
						pathingB11="cd ../../B11/Sing/"
						pathingC12="cd ../../C12/Sing/"
						pathingCa40="cd ../../Ca40/Sing/"
						pathingCa48="cd ../../Ca48/Sing/"
						pathingFe54="cd ../../Fe54/Sing/"
						pathingAu197="cd ../../Au197/Sing/"
						pathingD2="cd ../../D2/Sing/"
				fi
		elif [ $pass == 'pass4' ]
			then
				if [ $sorc == 'coin' ]
					then
						pathingBe9="cd Mf/pass4/Be9/Coin/"
						pathingB10="cd ../../B10/Coin/"
						pathingB11="cd ../../B11/Coin/"
						pathingC12="cd ../../C12/Coin/"
						pathingCa40="cd ../../Ca40/Coin/"
						pathingCa48="cd ../../Ca48/Coin/"
						pathingFe54="cd ../../Fe54/Coin/"
						pathingAu197="cd ../../Au197/Coin/"
						pathingD2="cd ../../D2/Coin/"
				elif [ $sorc == 'sing' ]
					then
						pathingBe9="cd Mf/pass4/Be9/Sing/"
						pathingB10="cd ../../B10/Sing/"
						pathingB11="cd ../../B11/Sing/"
						pathingC12="cd ../../C12/Sing/"
						pathingCa40="cd ../../Ca40/Sing/"
						pathingCa48="cd ../../Ca48/Sing/"
						pathingFe54="cd ../../Fe54/Sing/"
						pathingAu197="cd ../../Au197/Sing/"
						pathingD2="cd ../../D2/Sing/"
				fi
		fi
fi



eval ${pathingBe9}
targetBe9="Be9"
eval ${location}
run_rootBe9="root -b -l -q \"${file}( \\\"${run}\\\", \\\"${targetBe9}\\\", \\\"${pass}\\\", \\\"${sorc}\\\")\""
eval ${run_rootBe9}

eval ${pathingB10}
targetB10="B10"
eval ${location}
run_rootB10="root -b -l -q \"${file}( \\\"${run}\\\", \\\"${targetB10}\\\", \\\"${pass}\\\", \\\"${sorc}\\\")\""
eval ${run_rootB10}

eval ${pathingB11}
targetB11="B11"
eval ${location}
run_rootB11="root -b -l -q \"${file}( \\\"${run}\\\", \\\"${targetB11}\\\", \\\"${pass}\\\", \\\"${sorc}\\\")\""
eval ${run_rootB11}

eval ${pathingC12}
targetC12="C12"
eval ${location}
run_rootC12="root -b -l -q \"${file}( \\\"${run}\\\", \\\"${targetC12}\\\", \\\"${pass}\\\", \\\"${sorc}\\\")\""
eval ${run_rootC12}

eval ${pathingCa40}
targetCa40="Ca40"
eval ${location}
run_rootCa40="root -b -l -q \"${file}( \\\"${run}\\\", \\\"${targetCa40}\\\", \\\"${pass}\\\", \\\"${sorc}\\\")\""
eval ${run_rootCa40}

eval ${pathingCa48}
targetCa48="Ca48"
eval ${location}
run_rootCa48="root -b -l -q \"${file}( \\\"${run}\\\", \\\"${targetCa48}\\\", \\\"${pass}\\\", \\\"${sorc}\\\")\""
eval ${run_rootCa48}

eval ${pathingFe54}
targetFe54="Fe54"
eval ${location}
run_rootFe54="root -b -l -q \"${file}( \\\"${run}\\\", \\\"${targetFe54}\\\", \\\"${pass}\\\", \\\"${sorc}\\\")\""
eval ${run_rootFe54}

if [ $pass == 'pass2' ] || [ $pass == 'pass3' ] || [ $pass == 'pass4' ]
	then
		eval ${pathingAu197}
		targetAu197="Au197"
		eval ${location}
		run_rootAu197="root -b -l -q \"${file}( \\\"${run}\\\", \\\"${targetAu197}\\\", \\\"${pass}\\\", \\\"${sorc}\\\")\""
		eval ${run_rootAu197}
fi

if [ $pass == 'pass1' ] || [ $pass == 'pass3' ] || [ $pass == 'pass4' ]
	then
		eval ${pathingD2}
		targetD2="LD2"
		eval ${location}
		run_rootD2="root -b -l -q \"${file}( \\\"${run}\\\", \\\"${targetD2}\\\", \\\"${pass}\\\", \\\"${sorc}\\\")\""
		eval ${run_rootD2}
fi