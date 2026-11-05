C     NASA CR-176883: standard atmosphere, subroutine ATMCON (shared by
C     the three data banks; transcribed from the 16-series listing).
      SUBROUTINE ATMCON(ALTUDE,TEMP,PRES,DEN,VIS,WA)
      ALT=ALTUDE*.3048
      SLOPE=-.0065
      TEMP=288.16+SLOPE*ALT
      EXP=2.718281828
      TSEA=288.16
      PSEA=101325.
      PRES1=22700.
      GRAV=9.82
      R=288.
      POWER=GRAV/(SLOPE*R)
      PRES=PSEA*(1./(TEMP/TSEA))**POWER
      IF(ALT.LE.11000.) GO TO 1
      PRES=PRES1*EXP**(-GRAV/(R*216.66)*(ALT-11000.))
      TEMP=216.66
    1 WA=SQRT(1.4*R*TEMP)/.3048
      TEMP=(TEMP-273.15)*9./5.+491.67
      VIS=.35*(TEMP/492.)**1.5*(690./(TEMP+198.))*1.E-6
      PRES=PRES*2116.2/PSEA
      DEN=PRES/(1718.*TEMP)
      ALTUDE=ALT/.3048
      RETURN
      END
