import PySimpleGUI as sg
import run_simulation
from os import startfile
from Forces import choose_force_function as force


def bestParameters(dict):
    sg.theme('DarkAmber')

    # headings = ['Parameter', 'Best value']
    # header = [[sg.Text('  ')] + [sg.Text(h, size=(14, 1)) for h in headings]]
    #
    # input_rows = [[sg.Input(size=(15, 1), pad=(0, 0)) for col in range(4)] for
    #               row in range(10)]
    #
    # layout = header + input_rows
    #
    # window = sg.Window('Best parameters table', layout)
    # event, values = window.read()


    my_text = ''
    for i in dict:
        my_text = my_text + i + ': %s' % dict[i] + '\n'

    # All the stuff inside the window.
    layout = [[sg.Text('Here are the best parameters generated by the Monte Carlo simulation:')],
              [sg.popup_scrolled(my_text)],
              [sg.Submit(button_text='OK'), sg.Exit(button_text='Exit simulation')]]

    # Create the Window
    window = sg.Window('Best Parameters', layout)

    event, values = window.read()
    if event == 'OK':
        exit()
    elif event == 'Exit simulation':
        main()

    window.close()
    return


def MonteCarlo():
    """
    The Monte Carlo screen. Will guide the user which parameters he
    should bring to the program in order to run Monte Carlo simulation.
    The user can choose to go back to the opening screen, exit the program
    or press 'Run simulation'.
    :return: The path to the video generated by the program
    """
    sg.theme('DarkAmber')

    # All the stuff inside the window.
    layout = [[sg.Text('In order to run Monte Carlo simulation, please enter the following information:')],
              [sg.Text('All next attributes are MANDATORY:')],
              [sg.Text('Upload beginning state frame'), sg.InputText(),sg.FileBrowse()],
              [sg.Text('Upload end state frame'), sg.InputText(),sg.FileBrowse()],
              [sg.Text('Output folder'), sg.InputText(), sg.FolderBrowse()],
              [sg.Text('Enter Num of MC iterations: (need to be an integer)'), sg.InputText()],
              [sg.Text('Enter BD time delta: (in seconds)'), sg.InputText()],
              [sg.Text('Choose ONE of the force functions below: ')],
              [sg.Listbox(values=force.ALL_FORCE_FUNCTIONS, size=(30, 6))],
              [sg.Submit(button_text='Run simulation'), sg.Button(button_text="Back"), sg.Exit(button_text='Exit')]]

    # Create the Window
    window = sg.Window('Monte Carlo', layout)

    event, values = window.read()
    window.close()
    if event == 'Exit':
        exit()
    elif event == 'Back':
        main()
    elif event == 'Run simulation':
        start_frame, end_frame, output_folder, iteration_number, bd_deltaT, forceFunc = \
            values[0], values[1], values[2], int(values[3]), float(values[4]), values[5][0]


        output_tuple = run_simulation.run_mc_simulation(
            frame1=start_frame, frame2=end_frame, mc_iterations=iteration_number,
            bd_deltaT=bd_deltaT, output=output_folder,
            force_function_name=forceFunc)

        mov_path = output_tuple[0]
        vector_dict = output_tuple[1]

        try:
            bestParameters(vector_dict)
        except Exception as ex:
            pass

        return mov_path


def BrownianDynamics():
    """
    The Brownian Dynamics screen. Will guide the user which parameters he
    should bring to the program in order to run Brownian Dynamics simulation.
    The user can choose to go back to the opening screen, exit the program
    or press 'Run simulation'.
    :return: The path to the video generated by the program
    """
    sg.theme('DarkAmber')

    # All the stuff inside the window.
    layout = [[sg.Text('In order to run Brownian Dynamics simulation, please enter the following information:')],
              [sg.Text('Provide start frame OR frame sizes:')],
              [sg.Text('Upload beginning state frame'), sg.InputText(), sg.FileBrowse()],
              [sg.Text('Enter x-axis frame size: '), sg.InputText()],
              [sg.Text('Enter y-axis frame size: '), sg.InputText()],
              [sg.Text('Next attributes are MANDATORY:')],
              [sg.Text('Output folder'), sg.InputText(), sg.FolderBrowse()],
              [sg.Text('Enter Num of BD iterations:'), sg.InputText()],
              [sg.Text('Enter BD time delta: (in seconds)'), sg.InputText()],
              [sg.Text('Choose ONE of the force functions below: ')],
              [sg.Listbox(values=force.ALL_FORCE_FUNCTIONS, size=(30, 6))],
              [sg.Submit(button_text='Run simulation'), sg.Button(button_text="Back"), sg.Exit(button_text='Exit')]]

    # Create the Window
    window = sg.Window('Brownian Dynamics', layout)

    event, values = window.read()
    window.close()
    if event == 'Exit':
        exit()
    elif event == 'Back':
        main()
    elif event == 'Run simulation':
        start_frame, size_x, size_y, output_folder, iteration_number, bd_deltaT, forceFunc = \
            values[0], values[1], values[2], values[3], int(values[4]), float(values[5]), values[6][0]

        if start_frame is not "":
            mov_path = run_simulation.run_bd_simulation_based_on_frame(
                frame=start_frame, bd_num_iterations=iteration_number,
                bd_dt=bd_deltaT, output=output_folder,
                force_function_name=forceFunc)

        else:
            mov_path = run_simulation.run_bd_simulation_generate_conf(
                max_x=int(size_x), max_y=int(size_y),
                bd_num_iterations=iteration_number,
                bd_dt=bd_deltaT, output=output_folder,
                force_function_name=forceFunc)

        return mov_path


def openingFrame():
    """
    The opening screen of the program. Will show two possibilities to choose
    from- Brownian Dynamics simulation and Monte Carlo simulation. else, user
    can choose to exit the program.
    :return: The button pressed by the user
    """
    sg.theme('DarkAmber')

    # All the stuff inside the window.
    layout = [[sg.Text('Welcome to immune synapse dynamics simulation tool')],
              [sg.Text('Provided to you by 3D Data Processing in Structural Biology hackathon group number 1')],
              [sg.Text('©All rights reserved to Dror Bar, Lior Fishman, Tal Nisan, Chana Goldstein, Shani Cheskis and Liron Gershuny')],
              [sg.Text('Start by choosing the simulation you would like us to run for you: ')],
              [sg.Button(button_text="Brownian Dynamics"), sg.Button(button_text="Monte Carlo"),
               sg.Exit(button_text='Exit')]]

    # Create the Window
    window = sg.Window('Synapse simulation', layout)

    event = window.Read()[0]
    window.Close()
    return event


def main():
    event = openingFrame()

    if event == 'Brownian Dynamics':
        mov_path = BrownianDynamics()
    elif event == 'Monte Carlo':
        mov_path = MonteCarlo()
    elif event == 'Exit':
        exit()

    startfile(mov_path)
    main()


if __name__ == '__main__':
    main()
