Documentation protocol
==========================

This tutorial intends to help to develop and improve the documentation of the platform. 

Requirements
~~~~~~~~~~~~~~

To follow the next steps, the packages ``sphinx`` and ``sphinx_rtd_theme`` are needed. They can be installed executing the next commands:

	``pip install sphinx``

	``pip install sphinx_rtd_theme``


Steps
~~~~~~~~

The steps to follow are:

#. Fork the repository of pele_platform. See `pele_platform <https://github.com/NostrumBioDiscovery/pele_platform>`_.
	* Then clone the forked repository executing on your terminal:

	``git clone https://github.com/**your-github-user**/pele_platform.git``
	
#. Create a new branch to develop the documentation:

	``git branch -b PELE-docs``

#. **Create the .rst file**. Save it on the corresponding folder in the path */pele_platform/build/build_docs/source*, depending on the documentation page being created or improved, and start adding the information. See `the reStructuredText Cheat Sheet <cheatsheet_rst.html>`_ for basic RST syntax.

#. Once the rst file is ready, compile it into an html.
	* To create the file you just need to enter the following command on the *pele_platform/docs/build_docs* directory: ``make html``
#. **After the html file is correctly created, visualize it through your web browser of preference.** For example, you can execute this command to visualize the file on firefox: ``firefox /build/html/rest_of_the_path`` 
#. After checking the html file, push the results to your forked repository, as follows:

		``git add source``
	
		``git commit -m '[description of the updates made]'``
	
		``git push``

#. **As soon as the repository is updated you can start a pull request to merge the branches so the final documentation is visible in your page.** 
