function d = sanitizeNumericalOffsetsInSample(s,off)
% Given in OFF an offset and in S a sample of some distribution that has
% that offset, returns the same sample without the data that are
% indistinguishable from the offset.

    d = s(s > off);

end