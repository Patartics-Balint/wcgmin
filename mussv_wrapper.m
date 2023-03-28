function [bnd, info] = mussv_wrapper(matin, blk, opt, fixedBlkIdx)
	if verLessThan('matlab', '9.2') % R2016b or older
		undoc_input = {[], [], [], fixedBlkIdx};
	else
		undoc_input = {fixedBlkIdx};
	end
	[bnd, info] = mussv(matin, blk, opt, undoc_input{:});
end