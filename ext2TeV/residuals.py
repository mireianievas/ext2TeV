
class Residuals(object):
	def __init__(self):
		# dict with log of values of X and Y-residues for different type of sources.
		self.xlogval = {}
		self.ylogval = {}
		self.ylogmod = {}
		self.ylogres = {}
		self.ylogerr = {}

	def reset_cls(self,cls):
		self.xlogval[cls] = {}
		self.ylogval[cls] = {}
		self.ylogmod[cls] = {}
		self.ylogres[cls] = {}
		self.ylogerr[cls] = {}

	def reset_cls_model(self,cls,model):
		if cls not in self.xlogval:
			self.reset_cls(cls)

		self.xlogval[cls][model] = []
		self.ylogval[cls][model] = []
		self.ylogmod[cls][model] = []
		self.ylogres[cls][model] = []
		self.ylogerr[cls][model] = []
	

	def add_value(self, point):
		cls = point['srccls']
		model = point['model']
		xlogval = point['xlogval']
		ylogval = point['ylogval']
		ylogerr = point['ylogerr']
		ylogmod = point['ylogmod']
		# ylogres = point['ylogres']

		if cls not in self.xlogval:
			self.reset_cls(cls)

		if model not in self.xlogval[cls]:
			self.reset_cls_model(cls,model)

		self.xlogval[cls][model].append(xlogval)
		self.ylogval[cls][model].append(ylogval)
		self.ylogmod[cls][model].append(ylogmod)
		self.ylogres[cls][model].append((ylogval -ylogmod ) /ylogerr)
