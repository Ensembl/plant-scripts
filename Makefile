test_demo_user_scripts:
	cd demo_user_scripts && perl demo_test.t

clean_demo_user_scripts:	
	cd demo_user_scripts && rm -f Brachypodium*gz && rm -f Compara*gz
