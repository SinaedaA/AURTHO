IFS=$'\n'
TFscanout=$1
for query in `cat $TFscanout | awk '{print $3}' | uniq` 				# query = 3rd field in TF_scan_$strain.out, and always take only one
	do 
	echo "## query " $query " ##"
		#do
		#domain1=`grep -P "\s"$query"\s" $TFscanout | head -n1 | awk '{print $1}'`
	#echo "DEBUG: current query = $query" >&2
	domain1=`grep -m1 -e "$query" $TFscanout | awk '{print $1}'`
	#echo "DEBUG: current domain1 = $domain1" >&2
	echo $domain1
		case $domain1 in
		MazE_antitoxin) echo $query "AbrB"
			 ;;
                Phage_AlpA) echo $query "AlpA"
                         ;;
                Arg_repressor) echo $query "ArgR"
                         ;;
                HTH_5) echo $query "ArsR"
                         ;;
                HTH_CodY) echo $query "CodY"
                         ;;
                ComK) echo $query "ComK"
                         ;;
                Crl) echo $query "Crl"
                         ;;
                Trns_repr_metal) echo $query "CsoR"
                         ;;
                CtsR) echo $query "CtsR"
                         ;;
                HTH_DeoR) echo $query "DeoR"
                         ;;
                Smart00529) echo $query "DtxR"
                         ;;
                FaeA) echo $query "FaeA"
                         ;;
                Fe_dep_repress) echo $query "Fe_dep_repress"
                         ;;
                FeoC) echo $query "FeoC"
                         ;;
                FlhC) echo $query "CodY"
                         ;;
                FlhD) echo $query "FlhD"
                         ;;
                FUR) echo $query "Fur"
                         ;;
                GntR) echo $query "GntR"
                         ;;
                GutM) echo $query "GutM"
                         ;;
                Histone_HNS) echo $query "Hns"
                         ;;
                HrcA) echo $query "HrcA"
                         ;;
                HlxR) echo $query "HlxR"
                         ;;
                KorB) echo $query "KorB"
                         ;;
                LacI) echo $query "LacI"
                         ;;
                Lsr2) echo $query "Lsr2"
                         ;;
                HTH_1|LysR_substrate) echo $query "LysR"
                         ;;
                MarR) echo $query "MarR"
                         ;;
                MetJ) echo $query "MetJ"
                         ;;
                Mor) echo $query "Mor"
                         ;;
                MtlR) echo $query "MtlR"
                         ;;
                PadR) echo $query "PadR"
                         ;;
                PRD) echo $query "Prd"
                         ;;
                PucR|HTH_30) echo $query "PucR"
                         ;;
                PuR_N) echo $query "PuR"
                         ;;
                ROK) echo $query "Rok"
                         ;;
                ROS_MUCR) echo $query "Ros_MucR"
                         ;;
                HTH_6) echo $query "RpiR"
                         ;;
                Rrf2) echo $query "Rrf2"
                         ;;
                RtcR) echo $query "RtcR"
                         ;;
                SfsA) echo $query "SfsA"
                         ;;
                TetR_N) echo $query "TetR"
                         ;;
                TrmB) echo $query "TrmB"
                         ;;
                Trp_repressor) echo $query "TrpR"
                         ;;
                Whib) echo $query "WhiB"
                         ;;
                AraC_N) echo $query "AraC"
                         ;;
                AsnC_trans_reg|Smart00344) echo $query "AsnC"
                         ;;
                Bac_DNA_binding|Smart00411) echo $query "Bhl"
                         ;;
               # HTH_IclR|IclR) echo $query "IclR"
					 IclR) echo $query "IclR"
                         ;;
                GerE) echo $query "GerE"
                         ;;
                MerR|MerR-DNA-bind) echo $query "MerR"
                         ;;
                HTH_Mga|Mga) echo $query "Mga"
                         ;;
                MerR|MerR-DNA-bind) echo $query "MerR"
                         ;;
                BetR)
	                domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
			echo $domain2
			if [ "$domain2" == "Response_reg" ]
				then echo $query "BetR"
			fi
                         ;;
                CitT)
                        domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
                        echo $domain2
                        if [ "$domain2" == "Response_reg" ]
                                then echo $query "CitT"
                        fi
                         ;;
                Bac_DnaA|Bac_DnaA_C)
                        domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
			if [ "$domain2" == "$domain1" ] ; then continue; fi
                        echo $domain2
                        if [[ "$domain2" == "Bac_DnaA_C" || "$domain2" == "Bac_DnaA" ]]
                                then echo $query "DnaA"
                        fi
                         ;;
                Peptidase_S24|LexA_DNA_bind)
                        domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
			if [ "$domain2" == "$domain1" ] ; then continue; fi
                        echo $domain2
                        if [[ "$domain2" == "LexA_DNA_bind" || "$domain2" == "Peptidase_S24" ]]
                                then echo $query "LexA"
                        fi
                         ;;
                BTAD)
                        domains=`echo $(grep -e "\s"$query"\s" $TFscanout | awk '{print $1}')` ; echo $domains | grep -q "Trans_reg_C"	
                        if [ $? -eq 0 ]
                                then echo $query "Sarp"
                        fi
                         ;;
                Trans_reg_C)
                        domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
                        echo $domain2
			case $domain2 in
				BTAD) echo $query "Sarp"
					;;
				Response_reg) echo $query "OmpR"
					;;
			esac
                         ;;
                YcbB)
                        domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
                        echo $domain2
                        if [ "$domain"2 == "Response_reg" ]
                                then echo $query "YcbB"
                        fi
                         ;;
                LytTR)
                        domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
			echo $domain2
                        if [ "$domain2" == "Response_reg" ]
                                then echo $query "LyTR"
                        fi
                         ;;
                HTH_8)
                        domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
                        echo $domain2
                        case $domain2 in
                                Response_reg) 
					domain3=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n3 | tail -n1 | awk '{print $1}'`
					echo $domain3
					case $domain3 in
						AAA) echo $query "NtrC"
                                        		;;
						*) echo $query "PrrA"
							;;
					esac
					;;
                                AAA)
                                        domain3=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n3 | tail -n1 | awk '{print $1}'`
                                        echo $domain3
                        		if [ "$domain2" == "Response_reg" ]
		                                then echo $query "NtrC"
                        		fi
                                        ;;
				*) echo $query "Fis"
					;;
                        esac   
                         ;;
                AAA)
                        domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
                        echo $domain2
                        if [[ "$domain2" == "HTH_8" || "$domain2" == "Response_reg" ]]
                                then domain3=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n3 | tail -n1 | awk '{print $1}'`
				if [ "$domain3" == "$domain2" ] ; then continue; fi
				echo $domain3
				if [[ "$domain3" == "Response_Reg" || "$domain3" == "HTH_8" ]]
					then echo $query "NtrC"
				fi
                        fi
                         ;;
		HTH_3|Smart00530)
			domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
			if [ "$domain2" == "$domain1" ] ; then continue; fi
			echo $domain2
			case $domain2 in
				SinI) echo $query "SinR"
					;;
				Smart00530|HTH_3)
					domain3=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n3 | tail -n1 | awk '{print $1}'`
					if [ "$domain3" == "SinI" ]
						then echo $query "SinR"
					else
						echo $query "Xre" 
					fi
					;;
				*) echo $query "Xre"
					;;
			esac
			 ;;
		SinI)
			domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
                        echo $domain2
                        if [[ "$domain2" == "HTH_3" || "$domain2" == "Smart00530" ]]
				then echo $query "SinR"
			fi
			;;
                Crp|Smart00419)
                        domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
			if [ "$domain2" == "$domain1" ] ; then continue; fi
                        echo $domain2
                        case $domain2 in
                                Sugar-bind) echo $query "SorC"
                                        ;;
                                Smart00419|Crp)
                                        domain3=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n3 | tail -n1 | awk '{print $1}'`
                                        if [ "$domain3" == "Sugar-bind" ]
                                                then echo $query "SorC"
                                        else
                                                echo $query "Crp"
                                        fi
                                        ;;
                                *) echo $query "Crp"
                                        ;;
                        esac
                         ;;
                Sugar-bind)
                        domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
                        echo $domain2
                        if [[ "$domain2" == "Crp" || "$domain2" == "Smart00419" ]]
                                then echo $query "SorC"
                        fi
                        ;;
		Smart00421)
                        domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
                        echo $domain2
		        case $domain2 in
                        	GerE) echo $query "LuxR"
                                	;;
                                Response_reg) echo $query "NarL"
                                        ;;
                        esac
			;;	
                SpoOA_C)
                        domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
                        echo $domain2
                        if [ "$domain2" == "Response_reg" ]
                                then echo $query "SpoOA"
                        fi
                         ;;
		HTH_AraC)
                        domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
                        echo $domain2
                        case $domain2 in
                                Response_reg) echo $query "YesN"
                                        ;;
				*) echo $query "AraC"
                        esac
			;;
		Sigma54_AID|Sigma54_CBD|Sigma54_DBD)
                        domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
			if [ "$domain2" == "$domain1" ] ; then continue; fi
                        echo $domain2
                        if [[ "$domain2" == "Sigma54_CBD" || "$domain2" == "Sigma54_DBD" || "$domain2" == "Sigma54_AID" ]]
                                then domain3=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n3 | tail -n1 | awk '{print $1}'`
				if [ "$domain3" == "$domain2" ] ; then continue; fi
                                echo $domain3
                                if [[ "$domain3" == "Sigma54_DBD" || "$domain2" == "Sigma54_AID" || "$domain3" == "Sigma54_CBD" ]]
                                        then echo $query "RpoN"
                                fi
                        fi
                         ;;
                Sigma70_r3)
                        domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
                        echo $domain2
                        if [[ "$domain2" == "Sigma70_r2" || "$domain2" == "Sigma70_r4" ]]
                                then domain3=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n3 | tail -n1 | awk '{print $1}'`
				if [ "$domain3" == "$domain2" ] ; then continue; fi
                                echo $domain3
                                if [[ "$domain3" == "Sigma70_r4" || "$domain3" == "Sigma70_r2" ]]
                                        then echo $query "RpoD"
                                fi
                        fi
                         ;;
                Sigma70_r2|Sigma70_r4)
                        domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
			if [ "$domain2" == "$domain1" ] ; then continue; fi
                        echo $domain2
			case $domain2 in
				Sigma70_r4|Sigma70_r2)
					domain3=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n3 | tail -n1 | awk '{print $1}'`
	                                echo $domain3
					if [ "$domain3" == "Sigma70_r3" ]
						then echo $query "RpoD"
					else echo $query "Ecf"
					fi
					;;
				Sigma70_r3)
					domain3=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n3 | tail -n1 | awk '{print $1}'`
                                        echo $domain3
	                                if [[ "$domain3" == "Sigma70_r4" || "$domain3" == "Sigma70_r2" ]]
						then echo $query "RpoD"
					fi
					;;
			esac
			;;
                Response_reg)
                        domain2=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
                        echo $domain2
                        case $domain2 in
                                BetR) echo $query "BetR"
					;;
                                CitT) echo $query "CitT"
                                        ;;
                                LytTR) echo $query "LytTR"
                                        ;;
                                Smart00421) echo $query "NarL"
                                        ;;
                                Trans_reg_C) echo $query "OmpR"
                                        ;;
                                HTH_8) echo $query "PrrA"
                                        ;;
                                SpoOA_C) echo $query "SpoOA"
                                        ;;
                                YcbB) echo $query "YcbB"
                                        ;;
                                HTH_AraC) echo $query "YesN"
                                        ;;
				HTH_8|AAA)
					domain3=`grep -m1 -e "\s"$query"\s" $TFscanout | head -n2 | tail -n1 | awk '{print $1}'`
                        		echo $domain3
					if [ "$domain3" == "$domain2" ] ; then continue; fi
					if [[ "$domain3" == "AAA" || "$domain3" == "HTH_8" ]]
						then echo $query "NtrC"
					fi
					;;
			esac
		esac	

	#done

done
