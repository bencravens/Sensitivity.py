all:
	@python3 toy_system.py
plot:
    @python3 toy_system.py
    @okular *.png
log:
    @rm output.txt 
    @python3 toy_system.py >> output.txt

