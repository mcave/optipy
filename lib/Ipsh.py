import IPython
if IPython.version_info[0]>=1:
        from IPython.terminal.embed import InteractiveShellEmbed
else:
        from IPython.frontend.terminal.embed import InteractiveShellEmbed
from IPython.config.loader import Config
from inspect import currentframe

# Configure the prompt so that I know I am in a nested (embedded) shell
cfg = Config()
prompt_config = cfg.PromptManager
prompt_config.in_template = 'N.In <\\#>: '
prompt_config.in2_template = '   .\\D.: '
prompt_config.out_template = 'N.Out<\\#>: '
cfg.TerminalInteractiveShell.colors = 'LightBG'

# Messages displayed when I drop into and exit the shell.
banner_msg = ("Hit Ctrl-D to exit interpreter and continue program.\n"
	      "Note that if you use %kill_embedded, you can fully deactivate\n"
              "This embedded instance so it will never turn on again")   
exit_msg = ''

# Put ipshell() anywhere in your code where you want it to open.
ipshell = InteractiveShellEmbed(config=cfg, banner1=banner_msg, exit_msg=exit_msg)

def ipsh():
    frame = currentframe().f_back
    msg = 'Stopped at {0.f_code.co_filename} and line {0.f_lineno}'.format(frame)
    ipshell(msg,stack_depth=2) # Go back one level!

ipshell2 = InteractiveShellEmbed(config=cfg, banner1=banner_msg, exit_msg=exit_msg)

def ipsh2():
    frame = currentframe().f_back
    msg = 'Stopped at {0.f_code.co_filename} and line {0.f_lineno}'.format(frame)
    ipshell2(msg,stack_depth=2) # Go back one level!

ipshell3 = InteractiveShellEmbed(config=cfg, banner1=banner_msg, exit_msg=exit_msg)

def ipsh3():
    frame = currentframe().f_back
    msg = 'Stopped at {0.f_code.co_filename} and line {0.f_lineno}'.format(frame)
    ipshell3(msg,stack_depth=2) # Go back one level!
