################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/AlgorithmComputations.c \
../src/FileOperations.c \
../src/GeneticAlgorithm.c \
../src/MatrixComputations.c \
../src/QR_SearchingAlgorithm.c 

OBJS += \
./src/AlgorithmComputations.o \
./src/FileOperations.o \
./src/GeneticAlgorithm.o \
./src/MatrixComputations.o \
./src/QR_SearchingAlgorithm.o 

C_DEPS += \
./src/AlgorithmComputations.d \
./src/FileOperations.d \
./src/GeneticAlgorithm.d \
./src/MatrixComputations.d \
./src/QR_SearchingAlgorithm.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C Compiler'
	gcc -I"C:\cygwin64\usr\local\include" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


