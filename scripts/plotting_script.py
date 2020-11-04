import subprocess

plots = open("plots.txt", "r")

plot = plots.readline().rstrip()
commands = []
while plot:
    commands.append("python ../../plot_clustering.py"
                    + " " + plot + ".txt" + " "
                    +  plot + " " + "100")
    plot = plots.readline().rstrip()



for command in commands:
    subprocess.run(command)
