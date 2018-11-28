FUNCTION Volt_Convert, DN_value

; Function to convert DN values to a voltage [V]

    factor = -20d / 2d^16			; 20 V per 16 ADC bits (20./2.^16)

	return, factor * DN_value

END





FUNCTION Thermistor_RtoC, R

; Function to convert input thermistor resistance [*] to temperature [C]

; From George

; Initialize

	Tk = 273.15d

	a = 1.129241D-3						; manufacturer's values for a 10K resistor

	b = 2.341077D-4

	c = 8.775468D-8



	degrees_C = 1d / (a + b*alog (R) + C*alog (R)^3) - Tk

		; is this the solution from the Steinhart & Hart Eqn ignoring T^5 term?

	return, degrees_C

END



FUNCTION minTherm_Convert, DN_value

;SIM Focal plane temperatures (miniature thermistors)



    ; these values are from the sybase database table TimSimThermistorCalibration

    rref = 23340

    vref = 2d

    gain = 7.98d



	volts = Volt_Convert (DN_value)

    RT = rref * (Gain * vref + volts)/(Gain * vref - volts)  ; resistance of thermistor

	return, Thermistor_RtoC (RT)



END





FUNCTION Thermistor_Convert, DN_value, sima=sima, simb=simb

; Function to convert DN values to a thermistor temperature [C]



    ; these values are from the sybase database table TimSimThermistorCalibration

    if keyword_set(simb) then begin

        rtop = 19600d

        vtop = 7.1651d

    endif else begin

        rtop = 19600d

        vtop = 7.17d

    endelse



	volts = Volt_Convert (DN_value)

	RT = rtop * volts / (vtop - volts)	; resistance of thermistor

	return, Thermistor_RtoC (RT)

END



