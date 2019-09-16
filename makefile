all:
	@rm output.txt	
	@echo "enter desired sea ice concentration"
	@python3 toy_system.py >> output.txt
	@okular *.png &
