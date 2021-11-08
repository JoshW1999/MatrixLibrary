class NonSquareMatrixException extends Exception{

	NonSquareMatrixException(String message){
		super(message);
	}
}

class SingularMatrixException extends Exception{

	SingularMatrixException(String message){
		super(message);
	}
}

class IncompatibleDimensionsException extends Exception{
	IncompatibleDimensionsException(String message){
		super(message);
	}
}
