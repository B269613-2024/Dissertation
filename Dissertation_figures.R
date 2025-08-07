library(ggplot2)

####################################
########     Eco3R     #############
####################################

## Making dataframes for Eco3R HADDOCK metrics 

# Eco3R HADDOCK score dataframe 

HADDOCK_score_df_Eco3R<- data.frame(
  aptamer_version=c("parent", "tr1", "tr2", "tr3","tr4", "tr5","tr6", "tr7"),
  HADDOCK_score=c(-120, -99, -117.3, -113, -111.9, -106.3, -106.2, -120.2),
  stringsAsFactors = FALSE
)


# Eco3R Electrostatic energy dataframe 

Elec_E_df_Eco3R<- data.frame(
  aptamer_version=c("parent", "tr1", "tr2", "tr3","tr4", "tr5","tr6", "tr7"),
  Electrostatic_energy= c(-305.5, -70.8, -127.6, -193, -63.8, -110.8, -89.6, -133.1),
  stringsAsFactors = FALSE
)





# Eco3R Van Der Waals energy dataframe 
VanDerVaals_Eco3R<-data.frame(
  aptamer_version=c("parent", "tr1", "tr2", "tr3","tr4", "tr5","tr6", "tr7"),
  VanDerVaals_Energy=c(-44.5, -61.6, -65.4, -65.4, -68.4, -66.2, -64.2, -62.5),
  stringsAsFactors = FALSE
)


# Eco3R desolvation energy dataframe 

desolvation_df_Eco3R<-data.frame(
  aptamer_version=c("parent", "tr1", "tr2", "tr3","tr4", "tr5","tr6", "tr7"),
  desolvation_energy=c(-14.5, -23.3, -26.7, -9.1, -30.8,-18.1, -24.2, -31),
  stringsAsFactors = FALSE
)


## Assembling the four data frames into one long table

score_long_Eco3R <- rbind(
  data.frame(metric = "HADDOCK score",          colour = "seagreen",
             aptamer_version = HADDOCK_score_df_Eco3R$aptamer_version,
             value           = HADDOCK_score_df_Eco3R$HADDOCK_score),
  data.frame(metric = "Electrostatic energy",   colour = "yellow",
             aptamer_version = Elec_E_df_Eco3R$aptamer_version,
             value           = Elec_E_df_Eco3R$Electrostatic_energy),
  data.frame(metric = "Van der Waals energy",   colour = "orange",
             aptamer_version = VanDerVaals_Eco3R$aptamer_version,
             value           = VanDerVaals_Eco3R$VanDerVaals_Energy),
  data.frame(metric = "Desolvation energy",     colour = "pink",
             aptamer_version = desolvation_df_Eco3R$aptamer_version,
             value           = desolvation_df_Eco3R$desolvation_energy)
)

# Flagging the parent bars for easy visualization 

score_long_Eco3R$is_parent <- score_long_Eco3R$aptamer_version == "parent"

#  Plotting the four HADDOCK metrics 
dev.new()
ggplot(score_long_Eco3R,
       aes(x   = value,
           y   = aptamer_version,
           fill = interaction(metric, is_parent))) +     #Colours by the metric and parent is coloured red 
  geom_col() +
  geom_text(aes(label = value),                         #Adding labels (what score each bar has)
            hjust = 0.3, 
            size = 4) +
  facet_wrap(~ metric, scales = "free_y") +             
  scale_fill_manual(values = c(
    "HADDOCK score.FALSE"          = "seagreen",
    "HADDOCK score.TRUE"           = "red",
    "Electrostatic energy.FALSE"   = "yellow",
    "Electrostatic energy.TRUE"    = "red",
    "Van der Waals energy.FALSE"   = "orange",
    "Van der Waals energy.TRUE"    = "red",
    "Desolvation energy.FALSE"     = "pink",
    "Desolvation energy.TRUE"      = "red"
  ),
  guide = "none"  #No legend needed
  ) +
  labs(title = "Eco3R HADDOCK metrics",
       x = "Aptamer version",
       y = "Value") +
  theme_minimal()

## Eco3R hbonds

hbonds_Eco3R<- data.frame(
  aptamer_version=c("parent", "tr1", "tr2", "tr3", "tr4", "tr5", "tr6", "tr7"),
  hbond_number=c(10,4,11, 7, 6, 9, 6, 5 ),
  stringsAsFactors = FALSE
)
# Flagging the parent aptamer for visualization

hbonds_Eco3R$is_parent <- hbonds_Eco3R$aptamer_version == "parent"

# Hydrogen bond plot (across aptamer versions)

ggplot(hbonds_Eco3R,
       aes(x = aptamer_version,
           y = hbond_number,
           fill = is_parent)) +
  geom_col() +
  geom_text(aes(label = hbond_number),   #Adding label for hydrogen bond number 
            vjust = -0.5, 
            size = 4) +  
  scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "purple")) +   #legend: parent/trimmed
  labs(title = "Number of h-bonds in ECA1 aptamer variants",
       x = "Aptamer version",
       y = "Number of h-bonds",
       fill = "Parent") +
  theme_minimal() +
  ylim(0, max(hbonds_Eco3R$hbond_number) + 2)



####################################
########     ECA1     ##############
####################################

## Making dataframes for ECA1 HADDOCK metrics 

# ECA1 HADDOCK score dataframe 

HADDOCK_score_df_ECA1<- data.frame(
  aptamer_version=c("parent", "tr1", "tr2"),
  HADDOCK_score=c(-113.8, -126.4, -83.8),
  stringsAsFactors = FALSE
)

# ECA1 electrostatic energy dataframe 

Elec_E_df_ECA1<- data.frame(
  aptamer_version=c("parent", "tr1", "tr2"),
  Electrostatic_energy= c(-101.3, -140.3, -83.1 ),
  stringsAsFactors = FALSE
)
# ECA1 Van Der Waals energy dataframe

VanDerVaals_ECA1<-data.frame(
  aptamer_version=c("parent", "tr1", "tr2"),
  VanDerVaals_Energy=c(-66,-65.3, -52 ),
  stringsAsFactors = FALSE
)


# ECA1 desolvation energy dataframe

desolvation_df_ECA1<-data.frame(
  aptamer_version=c("parent", "tr1", "tr2"),
  desolvation_energy=c(-27.7,-33.1,-15.2 ),
  stringsAsFactors = FALSE
)

## Assembling the four data frames into one long table

score_long_ECA1 <- rbind(
  data.frame(metric = "HADDOCK score",          colour = "seagreen",
             aptamer_version = HADDOCK_score_df_ECA1$aptamer_version,
             value           = HADDOCK_score_df_ECA1$HADDOCK_score),
  data.frame(metric = "Electrostatic energy",   colour = "yellow",
             aptamer_version = Elec_E_df_ECA1$aptamer_version,
             value           = Elec_E_df_ECA1$Electrostatic_energy),
  data.frame(metric = "Van der Waals energy",   colour = "orange",
             aptamer_version = VanDerVaals_ECA1$aptamer_version,
             value           = VanDerVaals_ECA1$VanDerVaals_Energy),
  data.frame(metric = "Desolvation energy",     colour = "pink",
             aptamer_version = desolvation_df_ECA1$aptamer_version,
             value           = desolvation_df_ECA1$desolvation_energy)
)

# Flagging the parent bars for easy visualization
score_long$is_parent <- score_long$aptamer_version == "parent"

## Plotting the four HADDOCK metrics for ECA1 aptamer versions

dev.new()
ggplot(score_long,
       aes(x   = value,
           y   = aptamer_version,
           fill = interaction(metric, is_parent))) +    #Colours by the metric and parent is coloured red  
  geom_col() +
  geom_text(aes(label = value),                       #Adding labels (what score each bar has)
            hjust = -0.1, 
            size = 4) +
  facet_wrap(~ metric, scales = "free_y") +               
  scale_fill_manual(values = c(
    "HADDOCK score.FALSE"          = "seagreen",        #Assigning colours to parent and different HADDOCK metrics
    "HADDOCK score.TRUE"           = "red",
    "Electrostatic energy.FALSE"   = "yellow",
    "Electrostatic energy.TRUE"    = "red",
    "Van der Waals energy.FALSE"   = "orange",
    "Van der Waals energy.TRUE"    = "red",
    "Desolvation energy.FALSE"     = "pink",
    "Desolvation energy.TRUE"      = "red"
  ),
  guide = "none"   # no legend  because colour cues are obvious
 
  ) +
  labs(title = " HADDOCK metrics",
       x = "Aptamer version",
       y = " Value") +
  theme_minimal()



## ECA1 hbonds

# ECA1 hydrogen bond dataframe (by aptamer version)

hbonds_ECA1<- data.frame(
  aptamer_version=c("parent", "tr1", "tr2"),
  hbond_number=c(9, 7, 6),
  stringsAsFactors = FALSE
)

# Flagging the parent aptamer for visualization

hbonds_ECA1$is_parent <- hbonds_ECA1$aptamer_version == "parent"

# Hydrogen bond plot (across aptamer versions) 
dev.new()
ggplot(hbonds_ECA1,
       aes(x = aptamer_version,
           y = hbond_number,
           fill = is_parent)) +
  geom_col() +
  geom_text(aes(label = hbond_number),      #Adding label for hydrogen bond number
            vjust = -0.5, 
            size = 4) +  
  scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "purple")) +    #legend: parent/trimmed
  labs(title = "Number of h-bonds in ECA1 aptamer variants",
       x = "Aptamer version",
       y = "Number of h-bonds",
       fill = "Parent") +
  theme_minimal() +
  ylim(0, max(hbonds_ECA1$hbond_number) + 2)

####################################
########     ECA2     ##############
####################################

## Making dataframes for ECA2 HADDOCK metrics 

# ECA2 HADDOCK score dataframe 

HADDOCK_score_df_ECA2<- data.frame(
  aptamer_version=c("parent", "tr1", "tr2", "tr3"),
  HADDOCK_score=c(-110.2, -106.1, -123.9,-133.6),
  stringsAsFactors = FALSE
)

# ECA2 electrostatic energy dataframe

Elec_E_df_ECA2<- data.frame(
  aptamer_version=c("parent", "tr1", "tr2","tr3"),
  Electrostatic_energy= c(-248, -224.9,-186.4, -288.9 ),
  stringsAsFactors = FALSE
)

# ECA2 Van Der Waals energy dataframe

VanDerVaals_ECA2<-data.frame(
  aptamer_version=c("parent", "tr1", "tr2", "tr3"),
  VanDerVaals_Energy=c(-48.4,-49.2,-76.6,-77.7 ),
  stringsAsFactors = FALSE
)
# ECA2 desolvation energy dataframe

desolvation_df_ECA2<-data.frame(
  aptamer_version=c("parent", "tr1", "tr2", "tr3"),
  desolvation_energy=c(-12.3,-12,-10.1,1.7 ),
  stringsAsFactors = FALSE
)

## Assembling the four data frames into one long table

score_long_ECA2 <- rbind(
  data.frame(metric = "HADDOCK score",          colour = "seagreen",
             aptamer_version = HADDOCK_score_df_ECA2$aptamer_version,
             value           = HADDOCK_score_df_ECA2$HADDOCK_score),
  data.frame(metric = "Electrostatic energy",   colour = "yellow",
             aptamer_version = Elec_E_df_ECA2$aptamer_version,
             value           = Elec_E_df_ECA2$Electrostatic_energy),
  data.frame(metric = "Van der Waals energy",   colour = "orange",
             aptamer_version = VanDerVaals_ECA2$aptamer_version,
             value           = VanDerVaals_ECA2$VanDerVaals_Energy),
  data.frame(metric = "Desolvation energy",     colour = "pink",
             aptamer_version = desolvation_df_ECA2$aptamer_version,
             value           = desolvation_df_ECA2$desolvation_energy)
)

# Flagging the parent bars for easy visualization
score_long_ECA2$is_parent <- score_long_ECA2$aptamer_version == "parent"

## Plotting the four HADDOCK metrics for ECA2 aptamer versions

dev.new()
ggplot(score_long_ECA2,
       aes(x   = value,
           y   = aptamer_version,
           fill = interaction(metric, is_parent))) +      #Colours by the metric and parent is coloured red
  geom_col() +
  geom_text(aes(label = value),                           #Adding labels (what score each bar has)
            hjust = 0.3, 
            size = 4) +
  facet_wrap(~ metric, scales = "free_y") +               
  scale_fill_manual(values = c(
    "HADDOCK score.FALSE"          = "seagreen",          #Assigning colours to parent and different HADDOCK metrics
    "HADDOCK score.TRUE"           = "red",
    "Electrostatic energy.FALSE"   = "yellow",
    "Electrostatic energy.TRUE"    = "red",
    "Van der Waals energy.FALSE"   = "orange",
    "Van der Waals energy.TRUE"    = "red",
    "Desolvation energy.FALSE"     = "pink",
    "Desolvation energy.TRUE"      = "red"
  ),
  guide = "none"   # no legend  because colour cues are obvious
  
  ) +
  labs(title = " ECA2 HADDOCK metrics",
       x = "Aptamer version",
       y = " Value") +
  theme_minimal()

## ECA2 hbonds

#ECA2 hydrogen bond dataframe (by aptamer version)
hbonds_ECA2<- data.frame(
  aptamer_version=c("parent", "tr1", "tr2", "tr3"),
  hbond_number=c(11,6,9,15),
  stringsAsFactors = FALSE
)
# Flagging the parent aptamer for visualization
hbonds_ECA2$is_parent <- hbonds_ECA2$aptamer_version == "parent"

# Hydrogen bond plot (across aptamer versions) 
ggplot(hbonds_ECA2,
       aes(x = aptamer_version,
           y = hbond_number,
           fill = is_parent)) +
  geom_col() +
  geom_text(aes(label = hbond_number),                          #Adding label for hydrogen bond number
            vjust = -0.5, 
            size = 4) +  
  scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "purple")) +   #legend: parent/trimmed
  labs(title = "Number of h-bonds in ECA2 aptamer variants",
       x = "Aptamer version",
       y = "Number of h-bonds",
       fill = "Parent") +
  theme_minimal() +
  ylim(0, max(hbonds_ECA2$hbond_number) + 2)


####################################
########     Eco4F     #############
####################################

## Making dataframes for Eco4F HADDOCK metrics 

# Eco4F HADDOCK score dataframe

HADDOCK_score_df_Eco4F<- data.frame(
  aptamer_version=c("parent", "tr1", "tr2", "tr3"),
  HADDOCK_score=c(-150.8,-105.6, -141.8, -108.8),
  stringsAsFactors = FALSE
)

# Eco4F electrostatic energy dataframe

Elec_E_df_Eco4F<- data.frame(
  aptamer_version=c("parent", "tr1", "tr2","tr3"),
  Electrostatic_energy= c(-300.6, -58.1, -302.8,-144.4 ),
  stringsAsFactors = FALSE
)


# Eco4F Van Der Waals energy dataframe
VanDerVaals_Eco4F<-data.frame(
  aptamer_version=c("parent", "tr1", "tr2", "tr3"),
  VanDerVaals_Energy=c(-73.9,-58.8, -66.9,-70.7 ),
  stringsAsFactors = FALSE
)

# Eco4F desolvation energy dataframe
desolvation_df_Eco4F<-data.frame(
  aptamer_version=c("parent", "tr1", "tr2", "tr3"),
  desolvation_energy=c(-16.8,-35.3,-14.3,-9.2),
  stringsAsFactors = FALSE
)

## Assembling the four data frames into one long table

score_long_Eco4F <- rbind(
  data.frame(metric = "HADDOCK score",          colour = "seagreen",
             aptamer_version = HADDOCK_score_df_Eco4F$aptamer_version,
             value           = HADDOCK_score_df_Eco4F$HADDOCK_score),
  data.frame(metric = "Electrostatic energy",   colour = "yellow",
             aptamer_version = Elec_E_df_Eco4F$aptamer_version,
             value           = Elec_E_df_Eco4F$Electrostatic_energy),
  data.frame(metric = "Van der Waals energy",   colour = "orange",
             aptamer_version = VanDerVaals_Eco4F$aptamer_version,
             value           = VanDerVaals_Eco4F$VanDerVaals_Energy),
  data.frame(metric = "Desolvation energy",     colour = "pink",
             aptamer_version = desolvation_df_Eco4F$aptamer_version,
             value           = desolvation_df_Eco4F$desolvation_energy)
)

# Flagging the parent bars for easy visualization
score_long_Eco4F$is_parent <- score_long_Eco4F$aptamer_version == "parent"

## Plotting the four HADDOCK metrics for Eco4F aptamer versions

dev.new()
ggplot(score_long_Eco4F,
       aes(x   = value,
           y   = aptamer_version,
           fill = interaction(metric, is_parent))) +     #Colours by the metric and parent is coloured red
  geom_col() +
  geom_text(aes(label = value),                          #Adding labels (what score each bar has)
            hjust = 0.25, 
            size = 4) +
  facet_wrap(~ metric, scales = "free_y") +               
  scale_fill_manual(values = c(
    "HADDOCK score.FALSE"          = "seagreen",         #Assigning colours to parent and different HADDOCK metrics
    "HADDOCK score.TRUE"           = "red",
    "Electrostatic energy.FALSE"   = "yellow",
    "Electrostatic energy.TRUE"    = "red",
    "Van der Waals energy.FALSE"   = "orange",
    "Van der Waals energy.TRUE"    = "red",
    "Desolvation energy.FALSE"     = "pink",
    "Desolvation energy.TRUE"      = "red"
  ),
  guide = "none"   
  
  ) +
  labs(title = " Eco4F HADDOCK metrics",
       x = "Aptamer version",
       y = " Value") +
  theme_minimal()

## Eco4F hbonds

# Eco4F hydrogen bond dataframe (by aptamer version)
hbonds_Eco4F<- data.frame(
  aptamer_version=c("parent", "tr1", "tr2", "tr3"),
  hbond_number=c(19, 2, 12,11 ),
  stringsAsFactors = FALSE
)
# Flagging the parent aptamer for visualization
hbonds_Eco4F$is_parent <- hbonds_Eco4F$aptamer_version == "parent"

## Hydrogen bond plot (across aptamer versions) 

ggplot(hbonds_Eco4F,
       aes(x = aptamer_version,
           y = hbond_number,
           fill = is_parent)) +
  geom_col() +
  geom_text(aes(label = hbond_number),     #Adding label for hydrogen bond number
            vjust = -0.5, 
            size = 4) +  
  scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "purple")) +    #legend: parent/trimmed
  labs(title = "Number of h-bonds in Eco4F aptamer variants",
       x = "Aptamer version",
       y = "Number of h-bonds",
       fill = "Parent") +
  theme_minimal() +
  ylim(0, max(hbonds_Eco4F$hbond_number) + 2)
