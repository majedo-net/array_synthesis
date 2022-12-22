

class ArrayElementException(Exception):
    '''
    Exception for array element related errors
    '''
    def __init__(self, message="Element Error") -> None:
        self.message = message
        super().__init__(self.message)
