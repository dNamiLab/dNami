def color(str,message='com'):
	if message == 'com':
		return '\033[1;40;94m'+str+'\033[0m'
	elif message == 'error':
		return '\033[1;40;41m'+str+'\033[0m'

def exception(str,message='com'):
	if message == 'com':
		print(color('[info]'+str,message='com'))
	elif message == 'error':
		print(color('[error] '+str,message='error'))
		import sys
		sys.exit()

				