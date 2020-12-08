class Settings(object):

    def __init__(self):

        self._params = {
            'image_exts'    : [['.tiff', '.tif'], 'list'],

            }

    def get(self, key):
        return self._params[key][0]
