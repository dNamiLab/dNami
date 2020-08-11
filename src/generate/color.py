def color(str,message='com'):
	if message == 'com':
		return '\033[1;40;94m'+str+'\033[0m'
	elif message == 'error':
		return '\033[1;40;41m'+str+'\033[0m'
		