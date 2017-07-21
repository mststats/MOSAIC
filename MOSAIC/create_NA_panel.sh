for c in {1..22}
do
  paste panels/Moroccangenofile.$c panels/Druzegenofile.$c panels/Jordaniangenofile.$c \
        panels/Palestiniangenofile.$c panels/Tunisiangenofile.$c panels/Egyptiangenofile.$c \
	panels/Mozabitegenofile.$c panels/Bedouingenofile.$c -d "" > panels/NorthAfricangenofile.$c
done
