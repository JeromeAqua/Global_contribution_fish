disp("\begin{sidewaystable}[]")
tab = zeros(7,15);

tab(:,2) = TABLE(:,2); %seq. respi
tab(:,5) = TABLE(:,5); %seq. fecal
tab(:,8) = TABLE(:,11); %seq. carc
tab(:,11) = TABLE(:,8); %seq. OL
tab(:,14) = tab(:,2)+tab(:,5)+tab(:,8)+tab(:,11); %seq. tot

tab(:,1) = Tab_Activeflux(:,3); %active flux respi
tab(:,4) = Tab_Activeflux(:,1) + Tab_Passiveflux(:,1); %active + passive flux fecal pellets
tab(:,7) = Tab_Activeflux(:,2) + Tab_Passiveflux(:,2); %active + passive flux deadfalls
tab(:,10) = Tab_Activeflux(:,4); %active flux OL
tab(:,13) = tab(:,1) + tab(:,4) + tab(:,7) + tab(:,10);

tab(:,3) = tab(:,2) ./ tab(:,1);
tab(:,6) = tab(:,5) ./ tab(:,4);
tab(:,9) = tab(:,8) ./ tab(:,7);
tab(:,12) = tab(:,11) ./ tab(:,10);
tab(:,15) = tab(:,14) ./ tab(:,13);

disp("\caption{Total export, corresponding sequestration, and sequestration time for the different pathways considered in the model. Respiration pathway corresponds to animal respiration, fecal pellets pathway to bacterial respiration due to fecal pellet degradation, carcasses to natural mortality and other losses to all other losses. The range provided is obtained with the different scenarios of the sensitivity analysis.}")
disp("\label{table_results}")
disp("\begin{tabu}{|l||l|l|l||l|l|l||l|l|l||l|l|l||l|l|l|}")
disp("\hline")
disp("\rowfont{\scriptsize}")
disp("\multicolumn{1}{|c||}{}  & \multicolumn{3}{c||}{Respiration pathway}  & \multicolumn{3}{c||}{Fecal pellets pathway} & \multicolumn{3}{c||}{Carcasses pathway} & \multicolumn{3}{c||}{Other losses}   & \multicolumn{3}{c|}{Total}  \\ \cline{2-16} ")
disp("\rowfont{\scriptsize}")
disp("Organism & \begin{tabular}[c]{@{}l@{}}Export\\ {[}PgC  yr$^{-1}${]}\end{tabular} & \begin{tabular}[c]{@{}l@{}}Sequestr.\\ {[}PgC{]}\end{tabular} & \begin{tabular}[c]{@{}l@{}}Sequestr.\\ time {[}yr{]}\end{tabular} & \begin{tabular}[c]{@{}l@{}}Export\\  {[}PgC yr$^{-1}${]}\end{tabular} & \begin{tabular}[c]{@{}l@{}}Sequestr.\\ {[}PgC{]}\end{tabular} & \begin{tabular}[c]{@{}l@{}}Sequestr.\\ time {[}yr{]}\end{tabular} & \begin{tabular}[c]{@{}l@{}}Export\\  {[}PgC yr$^{-1}${]}\end{tabular} & \begin{tabular}[c]{@{}l@{}}Sequestr.\\ {[}PgC{]}\end{tabular} & \begin{tabular}[c]{@{}l@{}}Sequestr.\\ time {[}yr{]}\end{tabular}& \begin{tabular}[c]{@{}l@{}}Export\\  {[}PgC yr$^{-1}${]}\end{tabular} & \begin{tabular}[c]{@{}l@{}}Sequestr.\\ {[}PgC{]}\end{tabular} & \begin{tabular}[c]{@{}l@{}}Sequestr.\\ time {[}yr{]}\end{tabular} & \begin{tabular}[c]{@{}l@{}}Export\\  {[}PgC yr$^{-1}${]}\end{tabular} & \begin{tabular}[c]{@{}l@{}}Sequestr.\\ {[}PgC{]}\end{tabular} & \begin{tabular}[c]{@{}l@{}}Sequestr.\\ time {[}yr{]}\end{tabular}\\ \hline")
disp("\rowfont{\scriptsize}")
disp(strcat("\begin{tabular}[c]{@{}l@{}}Meso\\ zoopl.\end{tabular} & ", num2str(tab(1,1), '%.1f')," & ", num2str(tab(1,2), '%.1f')," & ", num2str(tab(1,3), '%.0f'),"  & ", num2str(tab(1,4), '%.1f'),"& ", num2str(tab(1,5), '%.0f'), "& ", num2str(tab(1,6), '%.0f'), " & ", num2str(tab(1,7), '%.1f'), " & ", num2str(tab(1,8), '%.0f'), " & ", num2str(tab(1,9), '%.0f'), " & ", num2str(tab(1,10), '%.0e'), " & ", num2str(tab(1,11), '%.2f'), " & ", num2str(tab(1,12), '%.0f') , "& ", num2str(tab(1,13), '%.1f'), " & ", num2str(tab(1,14), '%.0f'), " & ", num2str(tab(1,15), '%.0f'), "\\"))
% disp("\rowfont{\tiny}")
% disp("& (0.07 - 0.3) & (0.2-0.5) & (2-3) & (1.0 - 5.1) & (76-489) & (39-170) & (0.4 - 1.0) & (72-262) & (103-364) & (0.02 - 0.3) & (0.1-0.8) & (0.2-4) & (1.5 - 6.7) & (148-752) & (98-112) \\")
% disp("\rowfont{\scriptsize}")
disp(strcat("\begin{tabular}[c]{@{}l@{}}Macro\\ zoopl.\end{tabular} & ", num2str(tab(2,1), '%.1f')," & ", num2str(tab(2,2), '%.0f')," & ", num2str(tab(2,3), '%.0f'),"  & ", num2str(tab(2,4), '%.1f'),"& ", num2str(tab(2,5), '%.0f'), "& ", num2str(tab(2,6), '%.0f'), " & ", num2str(tab(2,7), '%.1f'), " & ", num2str(tab(2,8), '%.0f'), " & ", num2str(tab(2,9), '%.0f'), " & ", num2str(tab(2,10), '%.1f'), " & ", num2str(tab(2,11), '%.1f'), " & ", num2str(tab(2,12), '%.0f') , "& ", num2str(tab(2,13), '%.1f'), " & ", num2str(tab(2,14), '%.0f'), " & ", num2str(tab(2,15), '%.0f'), "\\")) 
% disp("\rowfont{\tiny}")
% disp("& (0.2 - 0.7) & (12-46) & (49-61) & (1.2 - 4.0) & (190-843) & (100-298) & (0.2 - 0.7) & (62-257) & (292-522) & (0.2 - 1.4) & (0-48)* & (0-35)*  & (1.8 - 6.8) & (264-1194) & (143-175)\\")
% disp("\rowfont{\scriptsize}")
disp(strcat("\begin{tabular}[c]{@{}l@{}}Meso-\\ pelagic\end{tabular}& ", num2str(tab(3,1), '%.1f')," & ", num2str(tab(3,2), '%.0f')," & ", num2str(tab(3,3), '%.0f'),"  & ", num2str(tab(3,4), '%.1f'),"& ", num2str(tab(3,5), '%.0f'), "& ", num2str(tab(3,6), '%.0f'), " & ", num2str(tab(3,7), '%.0e'), " & ", num2str(tab(3,8), '%.0f'), " & ", num2str(tab(3,9), '%.0f'), " & ", num2str(tab(3,10), '%.1f'), " & ", num2str(tab(3,11), '%.0f'), " & ", num2str(tab(3,12), '%.0f') , "& ", num2str(tab(3,13), '%.1f'), " & ", num2str(tab(3,14), '%.0f'), " & ", num2str(tab(3,15), '%.0f'), "\\"))
% disp("\rowfont{\tiny}")
% disp("& (0.1 - 1.0) & (14-137) & (133-139) & (0.2 - 1.5) & (67-535) & (225-554) & (2e-3 - 0.02) & (1-11) & (508-716) & (0.07 - 0.3) & (8-41) & (104-134) & (0.4 - 2.8) & (90-724) & (225-259)\\ ")
% disp("\rowfont{\scriptsize}")
disp(strcat("\begin{tabular}[c]{@{}l@{}}Forage\\ fish\end{tabular}& ", num2str(tab(4,1), '%.0e')," & ", num2str(tab(4,2), '%.2f')," & ", num2str(tab(4,3), '%.0f'),"  & ", num2str(tab(4,4), '%.1f'),"& ", num2str(tab(4,5), '%.0f'), "& ", num2str(tab(4,6), '%.0f'), " & ", num2str(tab(4,7), '%.0e'), " & ", num2str(tab(4,8), '%.0f'), " & ", num2str(tab(4,9), '%.0f'), " & ", num2str(tab(4,10), '%.1f'), " & ", num2str(tab(4,11), '%.1f'), " & ", num2str(tab(4,12), '%.0f') , "& ", num2str(tab(4,13), '%.1f'), " & ", num2str(tab(4,14), '%.0f'), " & ", num2str(tab(4,15), '%.0f'), "\\"))
% disp("\rowfont{\tiny}")
% disp("& (0.01 - 0.04) & (0.06-0.3) & (3-7) & (0.1 - 0.4) & (42-162) & (196-511) & (5e-3 - 0.01) & (2-7) & (468-654) & (0.04 - 0.1) & (0.08-0.4) & (2-4) & (0.2 - 0.6) & (44-170) & (220-283)\\")
% disp("\rowfont{\scriptsize}")
disp(strcat("\begin{tabular}[c]{@{}l@{}}Large\\ pelagic\end{tabular}& ", num2str(tab(5,1), '%.0e')," & ", num2str(tab(5,2), '%.1f')," & ", num2str(tab(5,3), '%.0f'),"  & ", num2str(tab(5,4), '%.2f'),"& ", num2str(tab(5,5), '%.0f'), "& ", num2str(tab(5,6), '%.0f'), " & ", num2str(tab(5,7), '%.0e'), " & ", num2str(tab(5,8), '%.0f'), " & ", num2str(tab(5,9), '%.0f'), " & ", num2str(tab(5,10), '%.2f'), " & ", num2str(tab(5,11), '%.0f'), " & ", num2str(tab(5,12), '%.0f') , "& ", num2str(tab(5,13), '%.2f'), " & ", num2str(tab(5,14), '%.0f'), " & ", num2str(tab(5,15), '%.0f'), "\\"))
% disp("\rowfont{\tiny}")
% disp("& (1e-3 - 5e-3) & (0.2-0.7) & (138-156) & (0.03 - 0.2) & (20-116) & (467-684) & (9e-4 - 3e-3) & (0.6-2)  &(587-655) & (0.05 - 0.3) & (8-49) & (161-164) & (0.08 - 0.5) & (29-168) & (273-363) \\")
% disp("\rowfont{\scriptsize}")
disp(strcat("Jellyfish& ", num2str(tab(6,1), '%.0e')," & ", num2str(tab(6,2), '%.0f')," & ", num2str(tab(6,3), '%.0f'),"  & ", num2str(tab(6,4), '%.2f'),"& ", num2str(tab(6,5), '%.0f'), "& ", num2str(tab(6,6), '%.0f'), " & ", num2str(tab(6,7), '%.2f'), " & ", num2str(tab(6,8), '%.0f'), " & ", num2str(tab(6,9), '%.0f'), " & ", num2str(tab(6,10), '%.2f'), " & ", num2str(tab(6,11), '%.0f'), " & ", num2str(tab(6,12), '%.0f') , "& ", num2str(tab(6,13), '%.1f'), " & ", num2str(tab(6,14), '%.0f'), " & ", num2str(tab(6,15), '%.0f'), "\\"))
% disp("\rowfont{\tiny}")
% disp(" & (0.03 - 0.1) & (5-19) & (142-156) & (0.08 - 0.7) & (35-332) & (300-634) & (0.04 - 0.1) & (21-63) & (388-628) & (0.1 - 1.3) & (6-87) & (62-67) & (0.3 - 2.2) & (67-501) & (223-228)\\")
% disp("\rowfont{\scriptsize}")
disp(strcat("\textbf{Total}& \textbf{", num2str(tab(7,1), '%.1f'),"} & \textbf{", num2str(tab(7,2), '%.0f'),"} & \textbf{", num2str(tab(7,3), '%.0f'),"}  & \textbf{", num2str(tab(7,4), '%.1f'),"}& \textbf{", num2str(tab(7,5), '%.0f'), "}& \textbf{", num2str(tab(7,6), '%.0f'), "} & \textbf{", num2str(tab(7,7), '%.1f'), "} & \textbf{", num2str(tab(7,8), '%.0f'), "} & \textbf{", num2str(tab(7,9), '%.0f'), "} & \textbf{", num2str(tab(7,10), '%.1f'), "} & \textbf{", num2str(tab(7,11), '%.0f'), "} & \textbf{", num2str(tab(7,12), '%.0f') , "}& \textbf{", num2str(tab(7,13), '%.1f'), "} & \textbf{", num2str(tab(7,14), '%.0f'), "} & \textbf{", num2str(tab(7,15), '%.0f'), "} \\"))
% disp("\rowfont{\tiny}")
% disp("& \textbf{(0.6 - 1.9)} & \textbf{(49-179)} & \textbf{(61-100)} & \textbf{(2.8 - 11.5)} & \textbf{(509-2160)} & \textbf{(113-302)} & \textbf{(0.7 - 1.8)} & \textbf{(169-562)} & \textbf{(183-430)} & \textbf{(0.7 - 3.5)} & \textbf{(35-207)} & \textbf{(44-73)} & \textbf{(4.6 - 18.8)} & \textbf{(762-3017)} & \textbf{(159-161)}  \\")
disp("\hline")
disp("\multicolumn{16}{l}{\footnotesize *Here, the lowest sequestration value  obtained for other losses was negative (and was truncated to 0), meaning that some parameter settings (very high mesopelagic fish biomass) are not compatible} \\")
disp(" \multicolumn{13}{l}{\footnotesize with a viable macro-zooplankton population.}")
disp("\end{tabu}")
disp("\end{sidewaystable}")